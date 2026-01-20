#include "physics/rocket.hpp"
#include <cmath>
#include <algorithm>

static constexpr double G0 = 9.80665; // m/s^2
static constexpr double R_E_KM = 6371.0;

Rocket::Rocket(const std::vector<Stage>& stages): st(stages){
    if(!st.empty()){
        cur=0;
        fuel = st[0].fuel_kg;
        dry  = st[0].dry_kg;
        mass = dry + fuel;
    }
}

void Rocket::set_state(const RocketState& s){
    x=s.x; y=s.y; z=s.z;
    vx=s.vx; vy=s.vy; vz=s.vz;
    mass=s.mass;
}

RocketState Rocket::state() const{
    return {x,y,z,vx,vy,vz,mass};
}

double Rocket::mdot() const{
    if(cur>=st.size()) return 0.0;
    if(st[cur].isp_s<=0.0) return 0.0;
    return st[cur].thrust_n/(st[cur].isp_s*G0); // kg/s
}

void Rocket::step(double dt_s, double mu_km3_s2,
                  const double* target_xyz_km,
                  const double* target_vxyz_km_s,
                  double lead_tau_s){
    sep=false;
    powered=false;
    if(is_dead) return;

    // gravity (km/s^2)
    const double r2 = x*x + y*y + z*z;
    const double r  = std::sqrt(std::max(1e-12, r2));
    double ax = -mu_km3_s2 * x / (r*r*r);
    double ay = -mu_km3_s2 * y / (r*r*r);
    double az = -mu_km3_s2 * z / (r*r*r);

    if(cur < st.size() && fuel > 0.0 && mass > 1e-9){
        powered=true;

        // thrust dir
        double tx= x/r, ty= y/r, tz= z/r; // radial default
        if(target_xyz_km){
            double px = target_xyz_km[0];
            double py = target_xyz_km[1];
            double pz = target_xyz_km[2];
            if(target_vxyz_km_s && lead_tau_s > 0.0){
                px += target_vxyz_km_s[0] * lead_tau_s;
                py += target_vxyz_km_s[1] * lead_tau_s;
                pz += target_vxyz_km_s[2] * lead_tau_s;
            }
            double dx = px - x;
            double dy = py - y;
            double dz = pz - z;
            double d  = std::sqrt(std::max(1e-12, dx*dx + dy*dy + dz*dz));
            tx = dx/d; ty = dy/d; tz = dz/d;
        }

        // if below surface, force radial-out thrust
        const double alt_km = r - R_E_KM;
        if(alt_km < 0.0){ tx=x/r; ty=y/r; tz=z/r; }

        const double a_thrust_km_s2 = (st[cur].thrust_n / mass) / 1000.0;
        ax += a_thrust_km_s2 * tx;
        ay += a_thrust_km_s2 * ty;
        az += a_thrust_km_s2 * tz;

        // burn
        const double dm = mdot() * dt_s;
        const double burn = std::min(fuel, std::max(0.0, dm));
        fuel -= burn;
        mass = dry + fuel;

        if(fuel <= 0.0){
            sep=true;
            cur++;
            if(cur < st.size()){
                dry  = st[cur].dry_kg;
                fuel = st[cur].fuel_kg;
                mass = dry + fuel;
            }else{
                powered=false; // coast
            }
        }
    }

    // semi-implicit Euler
    vx += ax * dt_s;
    vy += ay * dt_s;
    vz += az * dt_s;

    x += vx * dt_s;
    y += vy * dt_s;
    z += vz * dt_s;

    // clamp to surface (simple ground collision)
    const double r2n = x*x + y*y + z*z;
    const double rn  = std::sqrt(std::max(1e-12, r2n));
    if(rn < R_E_KM){
        const double ux=x/rn, uy=y/rn, uz=z/rn;
        x = ux*R_E_KM; y=uy*R_E_KM; z=uz*R_E_KM;

        const double vr = vx*ux + vy*uy + vz*uz;
        if(vr < 0.0){
            vx -= vr*ux; vy -= vr*uy; vz -= vr*uz;
        }
    }

    if(!std::isfinite(x) || !std::isfinite(vx) || !std::isfinite(mass)) is_dead=true;
}

void Rocket::force_cutoff(){
    powered=false;
    // Mark as "no more stages" so mdot() -> 0 and stage_sep stays false going forward
    cur = st.size();
    fuel=0.0;
    dry=mass; // whatever is left
}
