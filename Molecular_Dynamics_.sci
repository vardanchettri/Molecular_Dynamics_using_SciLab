clear;

//--------------------------------------------------
// FUNCTION: lattice position generator
//--------------------------------------------------
function [a1, a2, a3]=lattice_pos(index, Lx, Ly, Lz, lattice_const)
   l = 1;
   for i = -Lx/2:Lx/2
       x = i * lattice_const;
       for j = -Ly/2:Ly/2
           y = j * lattice_const;
           for k = -Lz/2:Lz/2
               z = k * lattice_const;
               if (abs(x)<Lx/2) & (abs(y)<Ly/2) & (abs(z)<Lz/2) & (l<=index) then
                   a1 = x; a2 = y; a3 = z;
                   l = l + 1;
               end
           end
       end
   end
endfunction

//--------------------------------------------------
// PARAMETERS
//--------------------------------------------------
npart = 1000;
Lx = 10; Ly = 10; Lz = 10;
lattice_const = 1.05;

temp = 1.0;
sigma = 1.0;
epsilon = 1.0;
rc = 2.5 * sigma;
rc2 = rc^2;

//--------------------------------------------------
// ARRAYS
//--------------------------------------------------
x  = zeros(1,npart);
y  = zeros(1,npart);
z  = zeros(1,npart);

vx = zeros(1,npart);
vy = zeros(1,npart);
vz = zeros(1,npart);

fx = zeros(1,npart);
fy = zeros(1,npart);
fz = zeros(1,npart);

//--------------------------------------------------
// INITIAL POSITIONS & VELOCITIES
//--------------------------------------------------
for i = 1:npart
    [x(i),y(i),z(i)] = lattice_pos(i,Lx,Ly,Lz,lattice_const);
    vx(i) = rand() - 0.5;
    vy(i) = rand() - 0.5;
    vz(i) = rand() - 0.5;
end

vx = vx - mean(vx);
vy = vy - mean(vy);
vz = vz - mean(vz);

//--------------------------------------------------
// INITIAL FORCE COMPUTATION
//--------------------------------------------------
epot = 0.0;

for i = 1:npart-1
    for j = i+1:npart

        dx = x(i) - x(j);
        dy = y(i) - y(j);
        dz = z(i) - z(j);

        dx = dx - Lx * round(dx/Lx);
        dy = dy - Ly * round(dy/Ly);
        dz = dz - Lz * round(dz/Lz);

        r2 = dx^2 + dy^2 + dz^2;

        if r2 < rc2 then
            invr2  = 1/r2;
            invr6  = invr2^3;
            invr12 = invr6^2;

            epot = epot + 4*epsilon*(invr12 - invr6);

            ffac = 24*epsilon*invr2*(2*invr12 - invr6);

            fx(i) = fx(i) + ffac*dx;
            fy(i) = fy(i) + ffac*dy;
            fz(i) = fz(i) + ffac*dz;

            fx(j) = fx(j) - ffac*dx;
            fy(j) = fy(j) - ffac*dy;
            fz(j) = fz(j) - ffac*dz;
        end
    end
end

//--------------------------------------------------
// TIME INTEGRATION
//--------------------------------------------------
dt = 0.005;
nsteps = 1000;

epot_arr = zeros(1,nsteps);
ekin_arr = zeros(1,nsteps);
etot_arr = zeros(1,nsteps);

fd = mopen("traj.xyz", "wt");

for step = 1:nsteps

    //---------------- POSITION UPDATE ----------------
    for i = 1:npart
        x(i) = x(i) + vx(i)*dt + 0.5*fx(i)*dt^2;
        y(i) = y(i) + vy(i)*dt + 0.5*fy(i)*dt^2;
        z(i) = z(i) + vz(i)*dt + 0.5*fz(i)*dt^2;

        x(i) = x(i) - Lx * round(x(i)/Lx);
        y(i) = y(i) - Ly * round(y(i)/Ly);
        z(i) = z(i) - Lz * round(z(i)/Lz);
    end

    //---------------- SAVE OLD FORCES ----------------
    fx_old = fx;
    fy_old = fy;
    fz_old = fz;

    //---------------- RESET FORCES ----------------
    fx(:) = 0; fy(:) = 0; fz(:) = 0;
    epot = 0.0;

    //---------------- FORCE RECOMPUTATION ----------------
    for i = 1:npart-1
        for j = i+1:npart

            dx = x(i) - x(j);
            dy = y(i) - y(j);
            dz = z(i) - z(j);

            dx = dx - Lx * round(dx/Lx);
            dy = dy - Ly * round(dy/Ly);
            dz = dz - Lz * round(dz/Lz);

            r2 = dx^2 + dy^2 + dz^2;

            if r2 < rc2 then
                invr2  = 1/r2;
                invr6  = invr2^3;
                invr12 = invr6^2;

                epot = epot + 4*epsilon*(invr12 - invr6);

                ffac = 24*epsilon*invr2*(2*invr12 - invr6);

                fx(i) = fx(i) + ffac*dx;
                fy(i) = fy(i) + ffac*dy;
                fz(i) = fz(i) + ffac*dz;

                fx(j) = fx(j) - ffac*dx;
                fy(j) = fy(j) - ffac*dy;
                fz(j) = fz(j) - ffac*dz;
            end
        end
    end

    //---------------- VELOCITY UPDATE ----------------
    for i = 1:npart
        vx(i) = vx(i) + 0.5*(fx_old(i) + fx(i))*dt;
        vy(i) = vy(i) + 0.5*(fy_old(i) + fy(i))*dt;
        vz(i) = vz(i) + 0.5*(fz_old(i) + fz(i))*dt;
    end

    //---------------- ENERGIES ----------------
    ekin = 0.0;
    for i = 1:npart
        ekin = ekin + 0.5*(vx(i)^2 + vy(i)^2 + vz(i)^2);
    end

    epot_arr(step) = epot;
    ekin_arr(step) = ekin;
    etot_arr(step) = epot + ekin;

    //---------------- XYZ OUTPUT ----------------
    mfprintf(fd, "%d\n", npart);
    mfprintf(fd, "Step %d\n", step);
    for i = 1:npart
        mfprintf(fd, "Ar  %f  %f  %f\n", x(i), y(i), z(i));
    end

end

mclose(fd);

//--------------------------------------------------
// ENERGY PLOT
//--------------------------------------------------
t = 1:nsteps;
clf();
plot(t, epot_arr, 'r');
plot(t, ekin_arr, 'b');
plot(t, etot_arr, 'k');
legend("Potential", "Kinetic", "Total");
xlabel("Time step");
ylabel("Energy");
title("Energy vs Time (Lennard-Jones MD)");
