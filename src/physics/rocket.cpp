#include "physics/rocket.hpp"
#include <cmath>
#include <algorithm>

static constexpr double G0      = 9.80665; // m/s^2
static constexpr double R_E_KM  = 6371.0;

Rocket::Rocket(const std::vector<Stage>& stages): st(stages){
    cur=0;
    sep=false;
    powered = !st.empty();
    if(!st.empty()){
        fuel = st[0].fuel_kg;
        dry  = st[0].dry_kg;
        mass = std::max(1.0, fuel + dry);
    }else{
        mass = 1.0;
    }
}

void Rocket::set_state(const RocketState& s){
    x=s.x; y=s.y; z=s.z;
    vx=s.vx; vy=s.vy; vz=s.vz;
    mass=std::max(1.0, s.mass);
}

RocketState Rocket::state() const{
    return {x,y,z,vx,vy,vz,mass};
}

double Rocket::mdot() const{
    if(cur>=st.size()) return 0.0;
    if(st[cur].isp_s<=0.0) return 0.0;
    return st[cur].thrust_n/(st[cur].isp_s*G0); // kg/s
}

void Rocket::step(double dt_s, double mu_km3_s2, const double* target_xyz_km){
    sep=false;

    // radius + unit vectors
    double r = std::sqrt(x*x+y*y+z*z);
    if(!(r>0.0)) r = R_E_KM;

    // hard clamp to surface + remove inward radial component (prevents "digging")
    if(r < R_E_KM){
        double ux=x/r, uy=y/r, uz=z/r;
        x = R_E_KM*ux; y = R_E_KM*uy; z = R_E_KM*uz;
        double vr = vx*ux + vy*uy + vz*uz;
        if(vr < 0.0){
            vx -= vr*ux; vy -= vr*uy; vz -= vr*uz;
        }
        r = R_E_KM;
    }

    double ux=x/r, uy=y/r, uz=z/r;

    // gravity accel (km/s^2)
    double aG = -mu_km3_s2/(r*r);
    double ax = aG*ux;
    double ay = aG*uy;
    double az = aG*uz;

    // thrust accel (km/s^2), with simple guidance:
    // - early: mostly radial-up
    // - later: blend toward target direction (if provided)
    double aTx=0.0,aTy=0.0,aTz=0.0;

    if(powered && cur < st.size() && fuel > 0.0 && st[cur].thrust_n > 0.0){
        double m = std::max(1.0, dry + fuel);

        // thrust magnitude: convert (m/s^2) -> (km/s^2)
        double aT = (st[cur].thrust_n / m) * 1e-3;

        // desired direction
        double dx=ux, dy=uy, dz=uz; // default: radial-up
        if(target_xyz_km){
            double tx = target_xyz_km[0] - x;
            double ty = target_xyz_km[1] - y;
            double tz = target_xyz_km[2] - z;
            double tn = std::sqrt(tx*tx+ty*ty+tz*tz);
            if(tn > 1e-9){
                tx/=tn; ty/=tn; tz/=tn;

                // ensure not pointing into Earth: if dot(dir,up) < 0, blend back to up
                double dotUp = tx*ux + ty*uy + tz*uz;
                if(dotUp < 0.0){
                    double k = std::min(1.0, -dotUp);
                    tx = (1.0-k)*tx + k*ux;
                    ty = (1.0-k)*ty + k*uy;
                    tz = (1.0-k)*tz + k*uz;
                    double tn2 = std::sqrt(tx*tx+ty*ty+tz*tz);
                    if(tn2>1e-9){ tx/=tn2; ty/=tn2; tz/=tn2; }
                }

                // time-based blend: first ~90s mostly up, then transition toward target
                // (dt-based since we may not have absolute t here; approximate by velocity magnitude)
                double vmag = std::sqrt(vx*vx+vy*vy+vz*vz);
                double blend = std::clamp((vmag - 0.5) / 6.0, 0.0, 1.0); // heuristic
                dx = (1.0-blend)*ux + blend*tx;
                dy = (1.0-blend)*uy + blend*ty;
                dz = (1.0-blend)*uz + blend*tz;
                double dn = std::sqrt(dx*dx+dy*dy+dz*dz);
                if(dn>1e-9){ dx/=dn; dy/=dn; dz/=dn; }
            }
        }

        aTx = aT*dx;
        aTy = aT*dy;
        aTz = aT*dz;

        // burn prop
        double dm = mdot()*dt_s;
        if(dm > fuel) dm = fuel;
        fuel -= dm;

        if(fuel <= 1e-9){
            sep=true;
            // stage sep: drop dry mass of spent stage, advance
            if(cur < st.size()){
                // drop current stage dry mass at sep
                // (keep mass from remaining stages by resetting dry/fuel below)
            }
            cur++;

            if(cur < st.size()){
                fuel = st[cur].fuel_kg;
                dry  = st[cur].dry_kg;
                powered=true;
            }else{
                powered=false; // coast, do NOT "die"
                fuel=0.0;
                dry=std::max(1.0, dry);
            }
        }

        mass = std::max(1.0, dry + fuel);
    }else{
        powered=false;
        mass = std::max(1.0, dry + fuel);
    }

    // integrate (semi-implicit Euler)
    ax += aTx; ay += aTy; az += aTz;

    vx += ax*dt_s;
    vy += ay*dt_s;
    vz += az*dt_s;

    x  += vx*dt_s;
    y  += vy*dt_s;
    z  += vz*dt_s;

    // clamp again after step
    r = std::sqrt(x*x+y*y+z*z);
    if(r < R_E_KM){
        double u2x=x/r, u2y=y/r, u2z=z/r;
        x = R_E_KM*u2x; y = R_E_KM*u2y; z = R_E_KM*u2z;
        double vr = vx*u2x + vy*u2y + vz*u2z;
        if(vr < 0.0){
            vx -= vr*u2x; vy -= vr*u2y; vz -= vr*u2z;
        }
    }
}
