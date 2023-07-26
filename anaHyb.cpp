// Created by Abhinav Singh - Copyright (c)
// ANASOL Interpolation case
#define SE_CLASS1
#include "config.h"
#define BOOST_MPL_CFG_NO_PREPROCESSED_HEADERS
#define BOOST_MPL_LIMIT_VECTOR_SIZE 40
#include <iostream>
#include "DCPSE/DCPSE_op/DCPSE_op.hpp"
#include "DCPSE/DCPSE_op/DCPSE_Solver.hpp"
#include "Operators/Vector/vector_dist_operators.hpp"
#include "Vector/vector_dist_subset.hpp"
#include "OdeIntegrators/OdeIntegrators.hpp"
#include <fstream>
#include "interpolation/mp4_kernel.hpp"
#include "interpolation/lambda_kernel.hpp"
#include "interpolation/z_spline.hpp"
#include "interpolation/interpolation.hpp"

void write_mag(std::string filename, double t,double mag){
    std::ofstream myFile(filename,std::ios_base::app);
    myFile <<t<<","<<mag<< "\n";
    myFile.close();
}

constexpr int x = 0;
constexpr int y = 1;
constexpr int z = 2;

double eta=1;
double nu=0;
double gama=1;
double Ks=1;
double Kb=1;
double Kt=1;
double lambda=-1;
double zeta=1;       //zeta=-1 corresponds to all previous notations before 2018 paper of Julicher. In our code zeta +ve would mean contractile and negative would mean extensile.
double timeTOL=1e-8;
double max_steady_tol=1e-8;
//dmu is set later as runtime parameter.
double dmu=1;
//Stress free.
double mu_t;
double mu_b;
double dt,V_err_eps;


int wr_f;
int wr_at;

constexpr int Polarization = 0;
constexpr int Velocity = 1;
constexpr int Vorticity = 2;
constexpr int ExtForce = 3;
constexpr int Pressure = 4;
constexpr int Strain_rate = 5;
constexpr int Stress = 6;
constexpr int MolField = 7;
constexpr int VRHS = 8;
constexpr int FE = 9;
constexpr int VT = 10;
constexpr int DV = 11;
constexpr int DPOL = 12;
constexpr int POLD = 13;
constexpr int DIV = 14;
constexpr int DELMU = 15;
constexpr int HPERP = 16;
constexpr int HPARR = 17;
constexpr int POLMAG = 18;
constexpr int TPOLMAG = 19;
constexpr int DT = 20;
constexpr int PARTICLENUMBER = 21;



void *vectorGlobal,*vectorGlobal_bulk,*vectorGlobal_boundary_up,*vectorGlobal_boundary_down,*global_gd,*global_vd;
const openfpm::vector<std::string> PropNAMES={"00-Polarization","01-Velocity","02-Vorticity","03-ExternalForce","04-Pressure","05-StrainRate","06-Stress","07-MolecularField","08-VelocityRHS","09-FranckEnergyDensity","10-V_t","11-dV","12-dPol","13-P_old","14-div","15-delmu","16-Hperp","17-Hparallel","18-PolMag","19-tPolMag","20-dt","21-particleNumber"};
typedef aggregate<double[3], VectorS<3, double>, double[3][3], VectorS<3, double>, double, double[3][3], double[3][3], VectorS<3, double>, VectorS<3, double>, double,VectorS<3, double>,VectorS<3, double>,VectorS<3, double>,VectorS<3, double>,double,double,VectorS<3, double>,double,double,double,double,int> ActiveGel3d;
typedef aggregate<double[3],double> interpolType;

typedef vector_dist_ws<3, double, ActiveGel3d> vector_type;
typedef vector_dist_subset<3, double, ActiveGel3d> vector_type2;
typedef grid_dist_id<3, double, interpolType> grid_type;
typedef vector_dist<3, double, interpolType> vd_type;

typedef interpolate<vd_type,grid_type,mp4_kernel<double>> inte_type2m;
typedef interpolate<vector_type,grid_type,z_kernel<double,1>> inte_type2p;      // first order z-spline results in a simple copy operation

// sets a particle property (vector) zero, used for remeshing
template<unsigned int prp> void set_zero_mygrid(grid_type & gr)
{
    auto it = gr.getDomainGhostIterator();

    // calculate the velocity from the curl of phi
    while (it.isNext())
    {
        auto key = it.get();

        gr.template get<prp>(key)[0] = 0.0;
        gr.template get<prp>(key)[1] = 0.0;
        gr.template get<prp>(key)[2] = 0.0;

        ++it;
    }
}

// sets a particle property (scalar) zero, used for remeshing
template<unsigned int prp> void set_zero_mygrid2(grid_type & gr)
{
    auto it = gr.getDomainGhostIterator();

    // calculate the velocity from the curl of phi
    while (it.isNext())
    {
        auto key = it.get();

        gr.template get<prp>(key) = 0.0;
        ++it;
    }
}

// this function calculates the Lagrangian derivative Dp/Dt given the velocity field from the CalcVelocity function
template<typename DX,typename DY,typename DZ,typename DXX,typename DXY,typename DXZ,typename DYY,typename DYZ,typename DZZ>
struct PolarEv
{
    DX &Dx;
    DY &Dy;
    DZ &Dz;
    DXX &Dxx;
    DXY &Dxy;
    DXZ &Dxz;
    DYY &Dyy;
    DYZ &Dyz;
    DZZ &Dzz;
    //Constructor
    PolarEv(DX &Dx,DY &Dy,DZ &Dz,DXX &Dxx,DXY &Dxy,DXZ &Dxz,DYY &Dyy,DYZ &Dyz,DZZ &Dzz):Dx(Dx),Dy(Dy),Dz(Dz),Dxx(Dxx),Dxy(Dxy),Dxz(Dxz),Dyy(Dyy),Dyz(Dyz),Dzz(Dzz)
    {}

    void operator()( const state_type_3d_ofp &X , state_type_3d_ofp &dxdt , const double t ) const
    {
        vector_type &Particles= *(vector_type *) vectorGlobal;
        vector_type2 &Particles_bulk= *(vector_type2 *) vectorGlobal_bulk;
        auto &v_cl = create_vcluster();
        auto Pol = getV<Polarization>(Particles);
        auto Pol_bulk = getV<Polarization>(Particles_bulk);
        auto V = getV<Velocity>(Particles);
        auto W = getV<Vorticity>(Particles);
        auto u = getV<Strain_rate>(Particles);
        auto h = getV<MolField>(Particles);
        auto dPol= getV<DPOL>(Particles);
        auto dPol_bulk= getV<DPOL>(Particles_bulk);
        auto delmu=getV<DELMU>(Particles);
        auto Hperp=getV<HPERP>(Particles);
        auto Hpar=getV<HPARR>(Particles);
        auto PolMag=getV<POLMAG>(Particles);
        auto tPolMag=getV<TPOLMAG>(Particles);
        Pol_bulk[x]=X.data.get<0>();
        Pol_bulk[y]=X.data.get<1>();
        Pol_bulk[z]=X.data.get<2>();
        Particles.ghost_get<Polarization>(SKIP_LABELLING);

        auto px=Pol[x];
        auto py=Pol[y];
        auto pz=Pol[z];

        texp_v<double> dxpx=Dx(Pol[x]),dxpy=Dx(Pol[y]),dxpz=Dx(Pol[z]),dypx=Dy(Pol[x]),dypy=Dy(Pol[y]),dypz=Dy(Pol[z]),dzpx=Dz(Pol[x]),dzpy=Dz(Pol[y]),dzpz=Dz(Pol[z]),
        dxxpx=Dxx(Pol[x]),dxxpy=Dxx(Pol[y]),dxxpz=Dxx(Pol[z]),dyypx=Dyy(Pol[x]),dyypy=Dyy(Pol[y]),dyypz=Dyy(Pol[z]),dzzpx=Dzz(Pol[x]),dzzpy=Dzz(Pol[y]),dzzpz=Dzz(Pol[z]),dxypx=Dxy(Pol[x]),
        dxypy=Dxy(Pol[y]),dxypz=Dxy(Pol[z]),dxzpx=Dxz(Pol[x]),dxzpy=Dxz(Pol[y]),dxzpz=Dxz(Pol[z]),dyzpx=Dyz(Pol[x]),dyzpy=Dyz(Pol[y]),dyzpz=Dyz(Pol[z]),
        polmaginv;

        h[x]=Ks*(dxxpx + dxypy + dxzpz) +
             Kb*((-dxypy - dxzpz + dyypx + dzzpx)*px*px + (-dxypy + dyypx)*py*py + (dypy*dzpx + dxpy*(dypz - 2*dzpy) + dypx*dzpy + dxpz*(-dypy - 2*dzpz) + 2*dzpx*dzpz)*pz + (-dxzpz + dzzpx)*pz*pz +
                 py*(dypz*dzpx + dxpz*(-2*dypz + dzpy) + dxpy*(-2*dypy - dzpz) + dypx*(2*dypy + dzpz) + (-dxypz - dxzpy + 2*dyzpx)*pz) +
                 px*(-dxpy*dxpy - dxpz*dxpz + dypx*dypx + dypz*dypz + dzpx*dzpx - 2*dypz*dzpy + dzpy*dzpy + (-dyzpz + dzzpy)*py + (dyypz - dyzpy)*pz)) +
             Kt*((-dxzpz + dzzpx)*py*py + (dxpz*dypy -  dypy*dzpx + dypx*(2*dypz -  dzpy) + dxpy*(-3*dypz + 2*dzpy))*pz + (- dxypy +  dyypx)*pz*pz + py*(-dypz*dzpx + dxpz*(2*dypz - 3*dzpy) + 2*dzpx*dzpy +  dxpy*dzpz -   dypx*dzpz + ( dxypz +  dxzpy - 2*dyzpx)*pz) +
                 px*(-2*dypz*dypz + 4*dypz*dzpy - 2*dzpy*dzpy + ( dyzpz -  dzzpy)*py + (- dyypz + dyzpy)*pz));

        h[y]= Ks*(dxypx + dyypy + dyzpz) +
              Kb*((dxxpy - dxypx)*px*px + (dxxpy - dxypx - dyzpz + dzzpy)*py*py + (dxpz*dypx + dxpy*dzpx - 2*dypx*dzpx + dxpx*(-dypz + dzpy) - 2*dypz*dzpz + 2*dzpy*dzpz)*pz + (-dyzpz + dzzpy)*pz*pz +
                  py*(dxpy*dxpy + dxpz*dxpz - dypx*dypx - dypz*dypz - 2*dxpz*dzpx + dzpx*dzpx + dzpy*dzpy + (dxxpz - dxzpx)*pz) +
                  px*(dxpx*(2*dxpy - 2*dypx) + dypz*dzpx + dxpz*(-2*dypz + dzpy) + dxpy*dzpz - dypx*dzpz + (-dxzpz + dzzpx)*py + (-dxypz + 2*dxzpy - dyzpx)*pz)) +
              Kt*((-dyzpz + dzzpy)*px*px + (-3*dxpz*dypx + dxpy*(2*dxpz - dzpx) + 2*dypx*dzpx + dxpx*(dypz - dzpy))*pz + (dxxpy - dxypx)*pz*pz + py*(-2*dxpz*dxpz + 4*dxpz*dzpx - 2*dzpx*dzpx + (-dxxpz + dxzpx)*pz) +
                  px*(-3*dypz*dzpx + dxpz*(2*dypz - dzpy) + 2*dzpx*dzpy - dxpy*dzpz + dypx*dzpz + (dxzpz - dzzpx)*py + (dxypz - 2*dxzpy + dyzpx)*pz));

        h[z]=Ks*(dxzpx + dyzpy + dzzpz) +
             Kb*((dxxpz - dxzpx)*px*px + (dyypz - dyzpy)*py*py + (dxpy*dxpy + dxpz*dxpz - 2*dxpy*dypx + dypx*dypx + dypz*dypz - dzpx*dzpx - dzpy*dzpy)*pz + (dxxpz - dxzpx + dyypz - dyzpy)*pz*pz +
                 py*(dxpz*dypx + dxpy*dzpx - 2*dypx*dzpx + dypy*(2*dypz - 2*dzpy) + dxpx*(dypz - dzpy) + (dxxpy - dxypx)*pz) +
                 px*(dxpz*dypy + dxpx*(2*dxpz - 2*dzpx) - dypy*dzpx + dxpy*(dypz - 2*dzpy) + dypx*dzpy + (2*dxypz - dxzpy - dyzpx)*py + (-dxypy + dyypx)*pz))+
             Kt*((dyypz - dyzpy)*px*px + (dxxpz - dxzpx)*py*py + (-2*dxpy*dxpy + 4*dxpy*dypx - 2*dypx*dypx)*pz + py*(-dxpz*dypx + dxpy*(2*dxpz - 3*dzpx) + 2*dypx*dzpx + dxpx*(-dypz + dzpy) + (-dxxpy + dxypx)*pz) +
                 px*(-dxpz*dypy + dypy*dzpx + dypx*(2*dypz - 3*dzpy) + dxpy*(-dypz + 2*dzpy) + (-2*dxypz + dxzpy + dyzpx)*py + (dxypy - dyypx)*pz));

        //Now defined as p cross h
        Hperp[x]=h[z]*Pol[y]-h[y]*Pol[z];
        Hperp[y]=h[x]*Pol[z]-h[z]*Pol[x];
        Hperp[z]=h[y]*Pol[x]-h[x]*Pol[y];

        tPolMag=Pol[x]*Pol[x]+Pol[y]*Pol[y]+Pol[z]*Pol[z];

        auto it = Particles.getDomainIterator();
        while (it.isNext()) {
            auto p = it.get();
            Particles.getProp<TPOLMAG>(p) = (Particles.getProp<TPOLMAG>(p) == 0) ? 1 : Particles.getProp<TPOLMAG>(p);
            ++it;
        }
        Particles.ghost_get<TPOLMAG>(SKIP_LABELLING);
        polmaginv=1.0/tPolMag;
        Hpar=-gama*(lambda*delmu-nu*polmaginv*(Pol[x]*Pol[x]*u[x][x] + 2*Pol[x]*Pol[y]*u[x][y] + 2*Pol[x]*Pol[z]*u[x][z] +  Pol[y]*Pol[y]*u[y][y] + 2*Pol[y]*Pol[z]*u[y][z] +  Pol[z]*Pol[z]*u[z][z])); 

        dPol = 0.0;
        dPol_bulk[x] = (Hpar*Pol[x] - Hperp[z]*Pol[y] + Hperp[y]*Pol[z])/gama-nu*(Pol[x]*u[x][x]+Pol[y]*u[x][y]+Pol[z]*u[x][z]) + lambda*Pol[x]*delmu + (W[x][x]*Pol[x]+W[x][y]*Pol[y]+W[x][z]*Pol[z]);
        dPol_bulk[y] = (Hpar*Pol[y] - Hperp[x]*Pol[z] + Hperp[z]*Pol[x])/gama-nu*(Pol[x]*u[y][x]+Pol[y]*u[y][y]+Pol[z]*u[y][z]) + lambda*Pol[y]*delmu + (W[y][x]*Pol[x]+W[y][y]*Pol[y]+W[y][z]*Pol[z]);
        dPol_bulk[z] = (Hpar*Pol[z] - Hperp[y]*Pol[x] + Hperp[x]*Pol[y])/gama-nu*(Pol[x]*u[z][x]+Pol[y]*u[z][y]+Pol[z]*u[z][z]) + lambda*Pol[z]*delmu + (W[z][x]*Pol[x]+W[z][y]*Pol[y]+W[z][z]*Pol[z]);
        dPol=dPol/sqrt(PolMag);
        dxdt.data.get<0>()=dPol[x];
        dxdt.data.get<1>()=dPol[y];
        dxdt.data.get<2>()=dPol[z];
    }
};


// calculates the velocity field given the current polarization field, includs remeshing
template<typename DX,typename DY,typename DZ,typename DXX,typename DXY,typename DXZ,typename DYY,typename DYZ,typename DZZ,typename i2m,typename i2p>
struct CalcVelocity
{

    DX &Dx;
    DY &Dy;
    DZ &Dz;
    DXX &Dxx;
    DXY &Dxy;
    DXZ &Dxz;
    DYY &Dyy;
    DYZ &Dyz;
    DZZ &Dzz;
    i2m &inte2m;
    i2p &inte2p;
    //Constructor
    int ctr;
    // int ctr;
    double t_old;

    //Constructor
    CalcVelocity(DX &Dx,DY &Dy,DZ &Dz,DXX &Dxx,DXY &Dxy,DXZ &Dxz,DYY &Dyy,DYZ &Dyz,DZZ &Dzz,i2m &inte2m,i2p &inte2p):Dx(Dx),Dy(Dy),Dz(Dz),Dxx(Dxx),Dxy(Dxy),Dxz(Dxz),Dyy(Dyy),Dyz(Dyz),Dzz(Dzz),inte2m(inte2m),inte2p(inte2p)
    {
        ctr = 0;
        t_old = -dt;

    }

    void operator() (state_type_3d_ofp &state, double t)
    {

        timer tt;
        vector_type &Particles= *(vector_type *) vectorGlobal;
        vector_type2 &Particles_bulk= *(vector_type2 *) vectorGlobal_bulk;
        vector_type2 &Particles_boundary_up= *(vector_type2 *) vectorGlobal_boundary_up;
        vector_type2 &Particles_boundary_down= *(vector_type2 *) vectorGlobal_boundary_down;

        grid_type &gd= *(grid_type *) global_gd;
        vd_type &vd= *(vd_type *) global_vd;

        auto Pos = getV<PROP_POS>(Particles);
        auto Pol = getV<Polarization>(Particles);
        auto Pol_old = getV<POLD>(Particles);
        auto vdPol = getV<0>(vd);
        auto vdPolMag = getV<1>(vd);
        auto V = getV<Velocity>(Particles);
        auto PolMag=getV<POLMAG>(Particles);
        auto tPolMag=getV<TPOLMAG>(Particles);
        auto dt_track = getV<DT>(Particles);
        auto dPol= getV<DPOL>(Particles);

        auto Pol_bulk = getV<Polarization>(Particles_bulk);
        auto P_bulk = getV<Pressure>(Particles_bulk);
        auto RHS_bulk = getV<VRHS>(Particles_bulk);

        auto & bulk = Particles_bulk.getIds();
        auto & Boundary_up = Particles_boundary_up.getIds();
        auto & Boundary_down = Particles_boundary_down.getIds();

        auto Pol_up = getV<Polarization>(Particles_boundary_up);
        auto Pol_down = getV<Polarization>(Particles_boundary_down);

        auto &v_cl = create_vcluster();

        // if time is unequal to zero and if we are not in the first time step after loading from a checkpoint
        if (t != 0) {
            tt.start();
            Pol_bulk[x]=state.data.get<0>();
            Pol_bulk[y]=state.data.get<1>();
            Pol_bulk[z]=state.data.get<2>();

            PolMag=(Pol[x]*Pol[x]+Pol[y]*Pol[y]+Pol[z]*Pol[z]);
            // normalize polarization
//            Pol[z]=0;
            Pol[x]=Pol[x]/sqrt(PolMag);
            Pol[y]=Pol[y]/sqrt(PolMag);
            Pol[z]=Pol[z]/sqrt(PolMag);
            Particles.ghost_get<Polarization,POLMAG>(SKIP_LABELLING);
            //PosVD = Pos + (t-t_old)*V;
             // method of images 
            vd.clear();
            double spacing_np=gd.getSpacing()[1];
            auto it = Particles.getDomainIterator();
            while(it.isNext())
            { auto p=it.get();
                //std::cout<<"ChecK: "<<p.getKey()<<","<<Particles.size_local()<<std::endl;
                double xp = Particles.getPos(p)[0];
                double yp = Particles.getPos(p)[1];
                double zp = Particles.getPos(p)[2];
                // Create new particles in the whole domain
                vd.add();
                vd.getLastPos()[x]=xp+(t-t_old)*Particles.getProp<1>(p)[x];   // prop<1>: Velocity
                vd.getLastPos()[y]=yp+(t-t_old)*Particles.getProp<1>(p)[y];
                vd.getLastPos()[z]=zp+(t-t_old)*Particles.getProp<1>(p)[z];
                vd.getLastProp<0>()[x]=Particles.getProp<0>(p)[x];
                vd.getLastProp<0>()[y]=Particles.getProp<0>(p)[y];
                vd.getLastProp<0>()[z]=Particles.getProp<0>(p)[z];
                vd.getLastProp<1>()=Particles.getProp<POLMAG>(p);
                // if within 3 grid spacings (+ tolerance due to numerics), create mirror particles
                if (yp>10.0-(3*spacing_np+spacing_np/4.0) && yp<10.0-spacing_np/8.0)
                {
                    vd.add();
                    vd.getLastPos()[x]=xp+(t-t_old)*Particles.getProp<1>(p)[x];
                    vd.getLastPos()[y]=(yp+(t-t_old)*Particles.getProp<1>(p)[y])+2.0*(10.0-(yp+(t-t_old)*Particles.getProp<1>(p)[y]));
                    vd.getLastPos()[z]=zp+(t-t_old)*Particles.getProp<1>(p)[z];
                    vd.getLastProp<0>()[0]=-Particles.getProp<0>(p)[0];
                    vd.getLastProp<0>()[1]=2.0-Particles.getProp<0>(p)[1];
                    vd.getLastProp<0>()[2]=-Particles.getProp<0>(p)[2];
                    vd.getLastProp<1>()=2.0-Particles.getProp<POLMAG>(p);
                }
                else if (yp<(3*spacing_np+spacing_np/8.0)&& yp>0.0+spacing_np/8.0)
                {
                    vd.add();
                    vd.getLastPos()[x]=xp+(t-t_old)*Particles.getProp<1>(p)[x];
                    vd.getLastPos()[y]=(yp+(t-t_old)*Particles.getProp<1>(p)[y])-2.0*(yp+(t-t_old)*Particles.getProp<1>(p)[y]);
                    vd.getLastPos()[z]=zp+(t-t_old)*Particles.getProp<1>(p)[z];
                    vd.getLastProp<0>()[0]=2.0-Particles.getProp<0>(p)[0];
                    vd.getLastProp<0>()[1]=-Particles.getProp<0>(p)[1];
                    vd.getLastProp<0>()[2]=-Particles.getProp<0>(p)[2];
                    vd.getLastProp<1>()=2.0-Particles.getProp<POLMAG>(p);
                }
            ++it;
            }

            vd.map();
            vd.ghost_get<Polarization,1>();
            // set gd zero
            set_zero_mygrid<0>(gd);
            // interpolate with M'4 from vd to grid nodes (gd)
            inte2m.template p2m<0,0>(vd,gd);
            gd.template ghost_put<add_,0>();
            gd.ghost_get<0>();
            Pol[0]=0;
            Pol[1]=0;
            Pol[2]=0;
            // interpolate with 1st order z-spline from grid nodes (gd) to particles
            inte2p.template m2p<0,0>(gd,Particles);
            Pol_up[0]=0;
            Pol_up[1]=1;
            Pol_up[2]=0;
            Pol_down[0]=1;
            Pol_down[1]=0;
            Pol_down[2]=0;

            state.data.get<0>()=Pol[x];
            state.data.get<1>()=Pol[y];
            state.data.get<2>()=Pol[z];

            tt.stop();
            if (v_cl.rank() == 0) {
                std::cout << "interpolation took " << tt.getwct() << " seconds." << std::endl;
                std::cout << "Time step " << ctr << " : " << t << " over." << std::endl;
                std::cout << "----------------------------------------------------------" << std::endl;
            }
            ctr++;
        }

        auto Dyx = Dxy;
        auto Dzy = Dyz;
        auto Dzx = Dxz;

        for (int j = 0; j < bulk.size(); j++) {
            auto p = bulk.get<0>(j);
            Particles.getProp<POLMAG>(p) = (Particles.getProp<POLMAG>(p) == 0) ? 1 : Particles.getProp<POLMAG>(p);
        }

        // pmag should be close to numerical precision due to interpolation and normalizaiton
        double pmag=0;
        for (int j = 0; j < bulk.size(); j++) {
            auto p = bulk.get<0>(j);
            if(fabs(1.0-Particles.getProp<POLMAG>(p))>pmag)
            {
                pmag=fabs(1-Particles.getProp<POLMAG>(p));
            }
        }

        v_cl.max(pmag);
        v_cl.execute();
        if(v_cl.rank()==0)
        {std::cout<<"Max Polarity Magnitude Deviation: "<<pmag<<std::endl;
        write_mag("PolMag",t,pmag);
        }
       
        auto W = getV<Vorticity>(Particles);
        auto g = getV<ExtForce>(Particles);
        auto P = getV<Pressure>(Particles);
        auto u = getV<Strain_rate>(Particles);
        auto sigma = getV<Stress>(Particles);
        auto h = getV<MolField>(Particles);
        auto RHS = getV<VRHS>(Particles);
        auto FranckEnergyDensity = getV<FE>(Particles);
        auto V_t = getV<VT>(Particles);
        auto dV  = getV<DV>(Particles);
        auto div =getV<DIV>(Particles);
        auto delmu=getV<DELMU>(Particles);
        auto Hperp=getV<HPERP>(Particles);
        auto Hpar=getV<HPARR>(Particles);
        auto g_up = getV<ExtForce>(Particles_boundary_up);
        auto g_down = getV<ExtForce>(Particles_boundary_down);
        auto V_up = getV<Velocity>(Particles_boundary_up);
        auto V_down = getV<Velocity>(Particles_boundary_down);


        // set up linear system solver
        petsc_solver<double> solverPetsc;
        solverPetsc.setSolver(KSPGMRES);
        solverPetsc.setRestart(5000);
        solverPetsc.setPreconditioner(PCNONE);

        auto px=Pol[x];
        auto py=Pol[y];
        auto pz=Pol[z];
        tt.start();
        Particles.ghost_get<Polarization>(SKIP_LABELLING);
        texp_v<double> dxpx=Dx(Pol[x]), dxpy=Dx(Pol[y]), dxpz=Dx(Pol[z]), dypx=Dy(Pol[x]), dypy=Dy(Pol[y]), dypz=Dy(Pol[z]),dzpx=Dz(Pol[x]),dzpy=Dz(Pol[y]),dzpz=Dz(Pol[z]),
                dxxpx=Dxx(Pol[x]),dxxpy=Dxx(Pol[y]),dxxpz=Dxx(Pol[z]),dyypx=Dyy(Pol[x]), dyypy=Dyy(Pol[y]), dyypz=Dyy(Pol[z]), dzzpx=Dzz(Pol[x]), dzzpy=Dzz(Pol[y]),
                dzzpz=Dzz(Pol[z]), dxypx=Dxy(Pol[x]), dxypy=Dxy(Pol[y]), dxypz=Dxy(Pol[z]), dxzpx=Dxz(Pol[x]), dxzpy=Dxz(Pol[y]), dxzpz=Dxz(Pol[z]), dyzpx=Dyz(Pol[x]), dyzpy=Dyz(Pol[y]),
                dyzpz=Dyz(Pol[z]), 
                //Note if mu is changing in space delmu has to go inside the derivative.
                dxmuqxx=delmu*Dx((Pol[x]*Pol[x]-1/3.0*(Pol[x]*Pol[x]+Pol[y]*Pol[y]+Pol[z]*Pol[z]))),
                dymuqxy=delmu*Dy(Pol[x]*Pol[y]) ,
                dzmuqxz=delmu*Dz(Pol[x]*Pol[z]) ,
                dxmuqyx=delmu*Dx(Pol[y]*Pol[x]),
                dymuqyy=delmu*Dy((Pol[y]*Pol[y]-1/3.0*(Pol[x]*Pol[x]+Pol[y]*Pol[y]+Pol[z]*Pol[z]))),
                dzmuqyz=delmu*Dz(Pol[y]*Pol[z]) ,
                dxmuqzx=delmu*Dx(Pol[z]*Pol[x]) ,
                dymuqzy=delmu*Dy(Pol[z]*Pol[y]),
                dzmuqzz=delmu*Dz((Pol[z]*Pol[z]-1/3.0*(Pol[x]*Pol[x]+Pol[y]*Pol[y]+Pol[z]*Pol[z]))),
                polmaginv,dPolMag;

        FranckEnergyDensity = 0.5*Ks*(dxpx + dypy + dzpz)*(dxpx + dypy + dzpz) +
                              0.5*Kt*((dypz - dzpy)*px + (-dxpz + dzpx)*py + (dxpy - dypx)*pz)*((dypz - dzpy)*px + (-dxpz + dzpx)*py + (dxpy - dypx)*pz) +
                              0.5*Kb*((-dxpz*px + dzpx*px - dypz*py + dzpy*py)*
                                      (-dxpz*px + dzpx*px - dypz*py + dzpy*py) +
                                      (dxpy*py - dypx*py + dxpz*pz - dzpx*pz)*
                                      (dxpy*py - dypx*py + dxpz*pz - dzpx*pz) +
                                      (-dxpy*px + dypx*px + dypz*pz - dzpy*pz)*
                                      (-dxpy*px + dypx*px + dypz*pz - dzpy*pz));

        h[x]=Ks*(dxxpx + dxypy + dxzpz) +
             Kb*((-dxypy - dxzpz + dyypx + dzzpx)*px*px + (-dxypy + dyypx)*py*py + (dypy*dzpx + dxpy*(dypz - 2*dzpy) + dypx*dzpy + dxpz*(-dypy - 2*dzpz) + 2*dzpx*dzpz)*pz + (-dxzpz + dzzpx)*pz*pz +
                 py*(dypz*dzpx + dxpz*(-2*dypz + dzpy) + dxpy*(-2*dypy - dzpz) + dypx*(2*dypy + dzpz) + (-dxypz - dxzpy + 2*dyzpx)*pz) +
                 px*(-dxpy*dxpy - dxpz*dxpz + dypx*dypx + dypz*dypz + dzpx*dzpx - 2*dypz*dzpy + dzpy*dzpy + (-dyzpz + dzzpy)*py + (dyypz - dyzpy)*pz)) +
             Kt*((-dxzpz + dzzpx)*py*py + (dxpz*dypy -  dypy*dzpx + dypx*(2*dypz -  dzpy) + dxpy*(-3*dypz + 2*dzpy))*pz + (- dxypy +  dyypx)*pz*pz + py*(-dypz*dzpx + dxpz*(2*dypz - 3*dzpy) + 2*dzpx*dzpy +  dxpy*dzpz -   dypx*dzpz + ( dxypz +  dxzpy - 2*dyzpx)*pz) +
                 px*(-2*dypz*dypz + 4*dypz*dzpy - 2*dzpy*dzpy + ( dyzpz -  dzzpy)*py + (- dyypz + dyzpy)*pz));

        h[y]= Ks*(dxypx + dyypy + dyzpz) +
              Kb*((dxxpy - dxypx)*px*px + (dxxpy - dxypx - dyzpz + dzzpy)*py*py + (dxpz*dypx + dxpy*dzpx - 2*dypx*dzpx + dxpx*(-dypz + dzpy) - 2*dypz*dzpz + 2*dzpy*dzpz)*pz + (-dyzpz + dzzpy)*pz*pz +
                  py*(dxpy*dxpy + dxpz*dxpz - dypx*dypx - dypz*dypz - 2*dxpz*dzpx + dzpx*dzpx + dzpy*dzpy + (dxxpz - dxzpx)*pz) +
                  px*(dxpx*(2*dxpy - 2*dypx) + dypz*dzpx + dxpz*(-2*dypz + dzpy) + dxpy*dzpz - dypx*dzpz + (-dxzpz + dzzpx)*py + (-dxypz + 2*dxzpy - dyzpx)*pz)) +
              Kt*((-dyzpz + dzzpy)*px*px + (-3*dxpz*dypx + dxpy*(2*dxpz - dzpx) + 2*dypx*dzpx + dxpx*(dypz - dzpy))*pz + (dxxpy - dxypx)*pz*pz + py*(-2*dxpz*dxpz + 4*dxpz*dzpx - 2*dzpx*dzpx + (-dxxpz + dxzpx)*pz) +
                  px*(-3*dypz*dzpx + dxpz*(2*dypz - dzpy) + 2*dzpx*dzpy - dxpy*dzpz + dypx*dzpz + (dxzpz - dzzpx)*py + (dxypz - 2*dxzpy + dyzpx)*pz));

        h[z]=Ks*(dxzpx + dyzpy + dzzpz) +
             Kb*((dxxpz - dxzpx)*px*px + (dyypz - dyzpy)*py*py + (dxpy*dxpy + dxpz*dxpz - 2*dxpy*dypx + dypx*dypx + dypz*dypz - dzpx*dzpx - dzpy*dzpy)*pz + (dxxpz - dxzpx + dyypz - dyzpy)*pz*pz +
                 py*(dxpz*dypx + dxpy*dzpx - 2*dypx*dzpx + dypy*(2*dypz - 2*dzpy) + dxpx*(dypz - dzpy) + (dxxpy - dxypx)*pz) +
                 px*(dxpz*dypy + dxpx*(2*dxpz - 2*dzpx) - dypy*dzpx + dxpy*(dypz - 2*dzpy) + dypx*dzpy + (2*dxypz - dxzpy - dyzpx)*py + (-dxypy + dyypx)*pz))+
             Kt*((dyypz - dyzpy)*px*px + (dxxpz - dxzpx)*py*py + (-2*dxpy*dxpy + 4*dxpy*dypx - 2*dypx*dypx)*pz + py*(-dxpz*dypx + dxpy*(2*dxpz - 3*dzpx) + 2*dypx*dzpx + dxpx*(-dypz + dzpy) + (-dxxpy + dxypx)*pz) +
                 px*(-dxpz*dypy + dypy*dzpx + dypx*(2*dypz - 3*dzpy) + dxpy*(-dypz + 2*dzpy) + (-2*dxypz + dxzpy + dyzpx)*py + (dxypy - dyypx)*pz));

        //Now defined as p cross h
        Hperp[x]=h[z]*Pol[y]-h[y]*Pol[z];
        Hperp[y]=h[x]*Pol[z]-h[z]*Pol[x];
        Hperp[z]=h[y]*Pol[x]-h[x]*Pol[y];

        sigma[x][x] =
                -dxpx*(dxpx + dypy + dzpz)*Ks -
                dxpy*Kt*pz*(dypz*px - dzpy*px - dxpz*py + dzpx*py + dxpy*pz - dypx*pz) -
                dxpz*Kt*py*(-dypz*px + dzpy*px + dxpz*py - dzpx*py - dxpy*pz + dypx*pz) -
                0.5*dxpz*Kb*(2*px*(dxpz*px - dzpx*px + (dypz - dzpy)*py) + 2*pz*(dxpy*py - dypx*py + (dxpz - dzpx)*pz)) -
                0.5*dxpy*Kb*(2*py*(dxpy*py - dypx*py + (dxpz - dzpx)*pz) + 2*px*(dxpy*px - dypx*px + (-dypz + dzpy)*pz));

        sigma[x][y] = 
            -Ks*dxpy*(dxpx + dypy + dzpz) + 
             Kt*(-px*dxpz + dxpx*pz)*(dypz*px - dzpy*px - dxpz*py + dzpx*py + dxpy*pz - dypx*pz) + 
             Kb*(-py*dxpz*(dxpz*px - dzpx*px + (dypz - dzpy)*py) + 
                dxpx*py*(dxpy*py - dypx*py + (dxpz - dzpx)*pz) - 
                dxpz*pz*(-px*dxpy + dypx*px + (dypz - dzpy)*pz) + 
                dxpx*px*(dxpy*px - dypx*px + (dzpy-dypz)*pz));

        sigma[x][z] =
                -dxpz*(dxpx + dypy + dzpz)*Ks +
                dxpy*Kt*px*(dypz*px - dzpy*px - dxpz*py + dzpx*py + dxpy*pz - dypx*pz) +
                dxpx*Kt*py*(-dypz*px + dzpy*px + dxpz*py - dzpx*py - dxpy*pz + dypx*pz) -
                0.5*dxpx*Kb*(2*px*(-dxpz*px + dzpx*px + (-dypz + dzpy)*py) - 2*pz*(dxpy*py - dypx*py + (dxpz - dzpx)*pz)) -
                0.5*dxpy*Kb*(2*py*(-dxpz*px + dzpx*px + (-dypz + dzpy)*py) - 2*pz*(-dxpy*px + dypx*px + (dypz - dzpy)*pz));
        sigma[y][x] =
                -dypx*(dxpx + dypy + dzpz)*Ks +
                dypz*Kt*py*(dypz*px - dzpy*px - dxpz*py + dzpx*py + dxpy*pz - dypx*pz) -
                dypy*Kt*pz*(dypz*px - dzpy*px - dxpz*py + dzpx*py + dxpy*pz - dypx*pz) -
                0.5*dypz*Kb*(2*px*(dxpz*px - dzpx*px + (dypz - dzpy)*py) + 2*pz*(dxpy*py - dypx*py + (dxpz - dzpx)*pz)) -
                0.5*dypy*Kb*(2*py*(dxpy*py - dypx*py + (dxpz - dzpx)*pz) + 2*px*(dxpy*px - dypx*px + (-dypz + dzpy)*pz));
        sigma[y][y] =
                -dypy*(dxpx + dypy + dzpz)*Ks -
                dypz*Kb*py*(dxpz*px - dzpx*px + (dypz - dzpy)*py) -
                dypz*Kt*px*(dypz*px - dzpy*px - dxpz*py + dzpx*py + dxpy*pz - dypx*pz) -
                dypx*Kt*pz*(-dypz*px + dzpy*px + dxpz*py - dzpx*py -dxpy*pz + dypx*pz) -
                dypx*Kb*py*(-dxpy*py + dypx*py + (-dxpz + dzpx)*pz) -
                dypx*Kb*px*(-dxpy*px + dypx*px + (dypz - dzpy)*pz) -
                dypz*Kb*pz*(-dxpy*px + dypx*px + (dypz - dzpy)*pz);
        sigma[y][z] =
                -dypz*(dxpx + dypy + dzpz)*Ks +
                dypy*Kt*px*(dypz*px - dzpy*px - dxpz*py + dzpx*py + dxpy*pz - dypx*pz) +
                dypx*Kt*py*(-dypz*px + dzpy*px + dxpz*py - dzpx*py -dxpy*pz + dypx*pz) -
                0.5*dypx*Kb*(2*px*(-dxpz*px + dzpx*px + (-dypz + dzpy)*py) - 2*pz*(dxpy*py - dypx*py + (dxpz - dzpx)*pz)) -
                0.5*dypy*Kb*(2*py*(-dxpz*px + dzpx*px + (-dypz + dzpy)*py) - 2*pz*(-dxpy*px + dypx*px + (dypz - dzpy)*pz));
        sigma[z][x] =
                -dzpx*(dxpx + dypy + dzpz)*Ks -
                dzpz*Kt*py*(-dypz*px + dzpy*px + dxpz*py - dzpx*py - dxpy*pz + dypx*pz) +
                dzpy*Kt*pz*(-dypz*px + dzpy*px + dxpz*py - dzpx*py - dxpy*pz + dypx*pz) -
                0.5*dzpz*Kb*(2*px*(dxpz*px - dzpx*px + (dypz - dzpy)*py) + 2*pz*(dxpy*py - dypx*py + (dxpz - dzpx)*pz)) -
                0.5*dzpy*Kb*(2*py*(dxpy*py - dypx*py + (dxpz - dzpx)*pz) + 2*px*(dxpy*px - dypx*px + (-dypz + dzpy)*pz));
        sigma[z][y] =
                -dzpy*(dxpx + dypy + dzpz)*Ks -
                dzpz*Kb*py*(dxpz*px - dzpx*px + (dypz - dzpy)*py) -
                dzpz*Kt*px*(dypz*px - dzpy*px - dxpz*py + dzpx*py + dxpy*pz - dypx*pz) +
                dzpx*Kt*pz*(dypz*px - dzpy*px - dxpz*py + dzpx*py + dxpy*pz - dypx*pz) -
                dzpx*Kb*py*(-dxpy*py + dypx*py + (-dxpz + dzpx)*pz) -
                dzpz*Kb*pz*(-dxpy*px + dypx*px + (dypz - dzpy)*pz) +
                dzpx*Kb*px*(dxpy*px - dypx*px + (-dypz + dzpy)*pz);

        sigma[z][z] =
                -dzpz*(dxpx + dypy + dzpz)*Ks -
                dzpx*Kt*py*(dypz*px - dzpy*px - dxpz*py + dzpx*py + dxpy*pz - dypx*pz) -
                dzpy*Kt*px*(-dypz*px + dzpy*px + dxpz*py - dzpx*py - dxpy*pz + dypx*pz) -
                0.5*dzpx*Kb*(2*px*(-dxpz*px + dzpx*px + (-dypz + dzpy)*py) - 2*pz*(dxpy*py - dypx*py + (dxpz - dzpx)*pz)) -
                0.5*dzpy*Kb*(2*py*(-dxpz*px + dzpx*px + (-dypz + dzpy)*py) - 2*pz*(-dxpy*px + dypx*px + (dypz - dzpy)*pz));

        Particles.ghost_get<Stress,MolField,HPERP>(SKIP_LABELLING);
        auto Hp1=Hperp[x];
        auto Hp2=Hperp[y];
        auto Hp3=Hperp[z];
        texp_v<double> 
        dxHp1=Dx(Hperp[x]),
        dxHp2=Dx(Hperp[y]),
        dxHp3=Dx(Hperp[z]),
        dyHp1=Dy(Hperp[x]),
        dyHp2=Dy(Hperp[y]),
        dyHp3=Dy(Hperp[z]),
        dzHp1=Dz(Hperp[x]),
        dzHp2=Dz(Hperp[y]),
        dzHp3=Dz(Hperp[z]);
        double dxmu=0;
        double dymu=0;
        double dzmu=0;

        // RHS (pressure is added later)
        dV[x]=  -(-0.5*dzpx*Hp2*px + 0.5*dypx*Hp3*px + 0.5*dzpx*Hp1*py + 0.5*dypy*Hp3*py + 0.5*dzpz*Hp3*py -
                 0.5*px*(-dzpy*Hp1 + dzpx*Hp2 + dzHp2*px - dzHp1*py) - 0.5*dypx*Hp1*pz - 0.5*dypy*Hp2*pz - 0.5*dzpz*Hp2*pz +
                0.5*px*(-dypz*Hp1 + dypx*Hp3 + dyHp3*px - dyHp1*pz) + 0.5*py*(-dypz*Hp2 + dypy*Hp3 + dyHp3*py - dyHp2*pz) -
                0.5*pz*(dzpz*Hp2 - dzpy*Hp3 - dzHp3*py + dzHp2*pz)
                +
                0.5*nu*(-(dypz*Hp1*px) + dzpy*Hp1*px + 2*dxpz*Hp2*px - 2*dzpx*Hp2*px - 2*dxpy*Hp3*px + 2*dypx*Hp3*px + dyHp3*px*px - dzHp2*px*px +
                    dzpx*Hp1*py + dypz*Hp2*py - 2*dxpx*Hp3*py - 2*dypy*Hp3*py - dzpz*Hp3*py - 2*dxHp3*px*py + dzHp1*px*py - dyHp3*py*py -
                    dypx*Hp1*pz + 2*dxpx*Hp2*pz + dypy*Hp2*pz + 2*dzpz*Hp2*pz - dzpy*Hp3*pz + 2*dxHp2*px*pz - dyHp1*px*pz + dyHp2*py*pz -
                    dzHp3*py*pz + dzHp2*pz*pz)
                + nu*gama*lambda*(-(dymu*px*py) - dzmu*px*pz + (dxmu*(-2*px*px + py*py + pz*pz))/3. +
                ((-4*dxpx*px)/3. - dypy*px - dzpz*px + (2*dxpy*py)/3. - dypx*py + (2*dxpz*pz)/3. - dzpx*pz)*delmu)
                + zeta*(dxmuqxx + dymuqxy + dzmuqxz) + Dx(sigma[x][x]) + Dy(sigma[x][y]) + Dz(sigma[x][z])
                );

        dV[y]=  -(-0.5*dzpy*Hp2*px - 0.5*dxpx*Hp3*px - 0.5*dzpz*Hp3*px + 0.5*dzpy*Hp1*py - 0.5*dxpy*Hp3*py +
                0.5*py*(dzpy*Hp1 - 1.*dzpx*Hp2 - 1.*dzHp2*px + dzHp1*py) + 0.5*dxpx*Hp1*pz + 0.5*dzpz*Hp1*pz + 0.5*dxpy*Hp2*pz -
                0.5*px*(-1.*dxpz*Hp1 + dxpx*Hp3 + dxHp3*px - 1.*dxHp1*pz) - 0.5*py*(-1.*dxpz*Hp2 + dxpy*Hp3 + dxHp3*py - 1.*dxHp2*pz) +
                0.5*pz*(dzpz*Hp1 - 1.*dzpx*Hp3 - 1.*dzHp3*px + dzHp1*pz)
                +
                0.5*nu*(-(dxpz*Hp1*px) - dzpy*Hp2*px + 2*dxpx*Hp3*px + 2*dypy*Hp3*px + dzpz*Hp3*px + dxHp3*px*px - 2*dypz*Hp1*py + 2*dzpy*Hp1*py +
                    dxpz*Hp2*py - dzpx*Hp2*py - 2*dxpy*Hp3*py + 2*dypx*Hp3*py + 2*dyHp3*px*py - dzHp2*px*py - dxHp3*py*py + dzHp1*py*py -
                    dxpx*Hp1*pz - 2*dypy*Hp1*pz - 2*dzpz*Hp1*pz + dxpy*Hp2*pz + dzpx*Hp3*pz - dxHp1*px*pz + dzHp3*px*pz + dxHp2*py*pz -
                    2*dyHp1*py*pz - dzHp1*pz*pz) +
                nu*gama*lambda*(-(dxmu*px*py) - dzmu*py*pz + (dymu*(px*px - 2*py*py + pz*pz))/3. +  (-(dxpy*px) + (2*dypx*px)/3. - dxpx*py - (4*dypy*py)/3. - dzpz*py + (2*dypz*pz)/3. - dzpy*pz)*delmu)
                +zeta*(dxmuqyx + dymuqyy + dzmuqyz) + Dx(sigma[y][x]) + Dy(sigma[y][y]) + Dz(sigma[y][z])
                )
                ;

        dV[z]= -(0.5*dxpx*Hp2*px + 0.5*dypy*Hp2*px + 0.5*dypz*Hp3*px - 0.5*dxpx*Hp1*py - 0.5*dypy*Hp1*py - 0.5*dxpz*Hp3*py +
                0.5*px*(-1.*dxpy*Hp1 + dxpx*Hp2 + dxHp2*px - 1.*dxHp1*py) - 0.5*py*(dypy*Hp1 - 1.*dypx*Hp2 - 1.*dyHp2*px + dyHp1*py) -
                0.5*dypz*Hp1*pz + 0.5*dxpz*Hp2*pz + 0.5*pz*(dxpz*Hp2 - 1.*dxpy*Hp3 - 1.*dxHp3*py + dxHp2*pz) -
                0.5*pz*(dypz*Hp1 - 1.*dypx*Hp3 - 1.*dyHp3*px + dyHp1*pz)
                +
                0.5*nu*(dxpy*Hp1*px - 2*dxpx*Hp2*px - dypy*Hp2*px - 2*dzpz*Hp2*px + dypz*Hp3*px - dxHp2*px*px + dxpx*Hp1*py + 2*dypy*Hp1*py +
                        2*dzpz*Hp1*py - dypx*Hp2*py - dxpz*Hp3*py + dxHp1*px*py - dyHp2*px*py + dyHp1*py*py - 2*dypz*Hp1*pz + 2*dzpy*Hp1*pz +
                        2*dxpz*Hp2*pz - 2*dzpx*Hp2*pz - dxpy*Hp3*pz + dypx*Hp3*pz + dyHp3*px*pz - 2*dzHp2*px*pz - dxHp3*py*pz + 2*dzHp1*py*pz +
                        dxHp2*pz*pz - dyHp1*pz*pz) +
                nu*gama*lambda*(-(dxmu*px*pz) - dymu*py*pz + (dzmu*(px*px + py*py - 2*pz*pz))/3. +
                                (-(dxpz*px) + (2*dzpx*px)/3. - dypz*py + (2*dzpy*py)/3. - dxpx*pz - dypy*pz - (4*dzpz*pz)/3.)*delmu)
                +zeta*(dxmuqzx + dymuqzy + dzmuqzz)+Dx(sigma[z][x])+Dy(sigma[z][y])+Dz(sigma[z][z])
                );

        Particles.ghost_get<DV>(SKIP_LABELLING);
        tt.stop();
        PolMag=Pol[x]*Pol[x]+Pol[y]*Pol[y]+Pol[z]*Pol[z];
        for (int j = 0; j < bulk.size(); j++) {
            auto p = bulk.get<0>(j);
            Particles.getProp<POLMAG>(p) = (Particles.getProp<POLMAG>(p) == 0) ? 1 : Particles.getProp<POLMAG>(p);
        }
        Particles.ghost_get<POLMAG>(SKIP_LABELLING);
        polmaginv=1.0/PolMag;

        auto zdxvx=Dx(V[x]);
        auto zdyvx=Dy(V[x]);
        auto zdzvx=Dz(V[x]);
        auto zdxxvx=Dxx(V[x]);
        auto zdxyvx=Dxy(V[x]);
        auto zdyxvx=Dxy(V[x]);
        auto zdxzvx=Dxz(V[x]);
        auto zdzxvx=Dxz(V[x]);
        auto zdyyvx=Dyy(V[x]);
        auto zdyzvx=Dyz(V[x]);
        auto zdzyvx=Dyz(V[x]);
        auto zdzzvx=Dzz(V[x]);
        auto zdxvy=Dx(V[y]);
        auto zdyvy=Dy(V[y]);
        auto zdzvy=Dz(V[y]);
        auto zdxxvy=Dxx(V[y]);
        auto zdxyvy=Dxy(V[y]);
        auto zdyxvy=Dxy(V[y]);
        auto zdxzvy=Dxz(V[y]);
        auto zdzxvy=Dxz(V[y]);
        auto zdyyvy=Dyy(V[y]);
        auto zdyzvy=Dyz(V[y]);
        auto zdzyvy=Dyz(V[y]);
        auto zdzzvy=Dzz(V[y]);
        auto zdxvz=Dx(V[z]);
        auto zdyvz=Dy(V[z]);
        auto zdzvz=Dz(V[z]);
        auto zdxxvz=Dxx(V[z]);
        auto zdxyvz=Dxy(V[z]);
        auto zdyxvz=Dxy(V[z]);
        auto zdxzvz=Dxz(V[z]);
        auto zdzxvz=Dxz(V[z]);
        auto zdyyvz=Dyy(V[z]);
        auto zdyzvz=Dyz(V[z]);
        auto zdzyvz=Dyz(V[z]);
        auto zdzzvz=Dzz(V[z]);

       auto Stokes1 = eta * (2*Dxx(V[x]) + Dxy(V[y]) + Dyy(V[x]) + Dxz(V[z]) + Dzz(V[x])) +
          (nu*nu*gama)*
          (
            polmaginv*polmaginv*(-2*px*py*(dypx*px + dypy*py + dypz*pz)*(px*px*zdxvx + px*py*(zdxvy + zdyvx) + py*py*zdyvy + px*pz*(zdxvz + zdzvx) +
                py*pz*(zdyvz + zdzvy) +pz*pz*zdzvz) -
             2*px*pz*(dzpx*px + dzpy*py + dzpz*pz)*(px*px*zdxvx + px*py*(zdxvy + zdyvx) + py*py*zdyvy + px*pz*(zdxvz + zdzvx) +
                py*pz*(zdyvz + zdzvy) +pz*pz*zdzvz) +
             1/3.0*(2*(dxpx*px + dxpy*py + dxpz*pz)*(-2*px*px + py*py + pz*pz)*
                (px*px*zdxvx + px*py*(zdxvy + zdyvx) + py*py*zdyvy + px*pz*(zdxvz + zdzvx) + py*pz*(zdyvz + zdzvy) +
                 pz*pz*zdzvz))) +
          polmaginv*(((4*dxpx*px)/3. + dypy*px + dzpz*px - (2*dxpy*py)/3. + dypx*py - (2*dxpz*pz)/3. + dzpx*pz)*
              (px*px*zdxvx + px*py*(zdxvy + zdyvx) + py*py*zdyvy + px*pz*(zdxvz + zdzvx) + py*pz*(zdyvz + zdzvy) +
               pz*pz*zdzvz) +
              1/3.0*((2*px*px - py*py - pz*pz)*
                (2*dxpx*px*zdxvx + px*px*zdxxvx + py*py*zdxxvy +pz*pz*zdxxvz + dxpy*px*(zdxvy + zdyvx) +
                  dxpx*py*(zdxvy + zdyvx) + 2*dxpy*py*zdyvy + px*py*(zdxxvy + zdyxvx) + dxpz*px*(zdxvz + zdzvx) +
                  dxpx*pz*(zdxvz + zdzvx) + dxpz*py*(zdyvz + zdzvy) + dxpy*pz*(zdyvz + zdzvy) + 2*dxpz*pz*zdzvz +
                  px*pz*(zdxxvz + zdzxvx) + py*pz*(zdyxvz + zdzxvy))) +
             px*py*(2*dypx*px*zdxvx + px*px*zdxyvx + dypy*px*(zdxvy + zdyvx) + dypx*py*(zdxvy + zdyvx) + 2*dypy*py*zdyvy +
                px*py*(zdxyvy + zdyyvx) + py*py*zdyyvy +pz*pz*zdyzvz + dypz*px*(zdxvz + zdzvx) + dypx*pz*(zdxvz + zdzvx) +
                dypz*py*(zdyvz + zdzvy) + dypy*pz*(zdyvz + zdzvy) + 2*dypz*pz*zdzvz + px*pz*(zdxyvz + zdzyvx) +
                py*pz*(zdyyvz + zdzyvy)) + px*pz*(2*dzpx*px*zdxvx + px*px*zdxzvx + dzpy*px*(zdxvy + zdyvx) +
                dzpx*py*(zdxvy + zdyvx) + 2*dzpy*py*zdyvy + px*py*(zdxzvy + zdyzvx) + py*py*zdyzvy +
                dzpz*px*(zdxvz + zdzvx) + dzpx*pz*(zdxvz + zdzvx) + dzpz*py*(zdyvz + zdzvy) + dzpy*pz*(zdyvz + zdzvy) +
                2*dzpz*pz*zdzvz + px*pz*(zdxzvz + zdzzvx) + py*pz*(zdyzvz + zdzzvy) +pz*pz*zdzzvz))
     );

        auto Stokes2 = eta * (2*Dyy(V[y]) + Dxx(V[y]) + Dxy(V[x]) + Dyz(V[z]) + Dzz(V[y]))+
        (gama*nu*nu)*(
            polmaginv*polmaginv*(-2*px*py*(dxpx*px + dxpy*py + dxpz*pz)*(px*px*zdxvx + px*py*(zdxvy + zdyvx) + py*py*zdyvy + px*pz*(zdxvz + zdzvx) +
                py*pz*(zdyvz + zdzvy) +pz*pz*zdzvz) -
             2*py*pz*(dzpx*px + dzpy*py + dzpz*pz)*(px*px*zdxvx + px*py*(zdxvy + zdyvx) + py*py*zdyvy + px*pz*(zdxvz + zdzvx) +
                py*pz*(zdyvz + zdzvy) +pz*pz*zdzvz) +
             1/3.0*(2*(dypx*px + dypy*py + dypz*pz)*(px*px - 2*py*py + pz*pz)*
                (px*px*zdxvx + px*py*(zdxvy + zdyvx) + py*py*zdyvy + px*pz*(zdxvz + zdzvx) + py*pz*(zdyvz + zdzvy) +
                 pz*pz*zdzvz)))+
          polmaginv*((dxpy*px - (2*dypx*px)/3. + dxpx*py + (4*dypy*py)/3. + dzpz*py - (2*dypz*pz)/3. + dzpy*pz)*
              (px*px*zdxvx + px*py*(zdxvy + zdyvx) + py*py*zdyvy + px*pz*(zdxvz + zdzvx) + py*pz*(zdyvz + zdzvy) +
               pz*pz*zdzvz) + px*py*(2*dxpx*px*zdxvx + px*px*zdxxvx + py*py*zdxxvy +pz*pz*zdxxvz + dxpy*px*(zdxvy + zdyvx) +
                dxpx*py*(zdxvy + zdyvx) + 2*dxpy*py*zdyvy + px*py*(zdxxvy + zdyxvx) + dxpz*px*(zdxvz + zdzvx) +
                dxpx*pz*(zdxvz + zdzvx) + dxpz*py*(zdyvz + zdzvy) + dxpy*pz*(zdyvz + zdzvy) + 2*dxpz*pz*zdzvz +
                px*pz*(zdxxvz + zdzxvx) + py*pz*(zdyxvz + zdzxvy)) +
             1/3.0*((-px*px + 2*py*py - pz*pz)*(2*dypx*px*zdxvx + px*px*zdxyvx + dypy*px*(zdxvy + zdyvx) + dypx*py*(zdxvy + zdyvx) +
                  2*dypy*py*zdyvy + px*py*(zdxyvy + zdyyvx) + py*py*zdyyvy +pz*pz*zdyzvz + dypz*px*(zdxvz + zdzvx) +
                  dypx*pz*(zdxvz + zdzvx) + dypz*py*(zdyvz + zdzvy) + dypy*pz*(zdyvz + zdzvy) + 2*dypz*pz*zdzvz +
                  px*pz*(zdxyvz + zdzyvx) + py*pz*(zdyyvz + zdzyvy))) +
             py*pz*(2*dzpx*px*zdxvx + px*px*zdxzvx + dzpy*px*(zdxvy + zdyvx) + dzpx*py*(zdxvy + zdyvx) + 2*dzpy*py*zdyvy +
                px*py*(zdxzvy + zdyzvx) + py*py*zdyzvy + dzpz*px*(zdxvz + zdzvx) + dzpx*pz*(zdxvz + zdzvx) +
                dzpz*py*(zdyvz + zdzvy) + dzpy*pz*(zdyvz + zdzvy) + 2*dzpz*pz*zdzvz + px*pz*(zdxzvz + zdzzvx) +
                py*pz*(zdyzvz + zdzzvy) +pz*pz*zdzzvz))
          );

        auto Stokes3 = eta * (2*Dzz(V[z]) + Dxx(V[z]) + Dxz(V[x]) + Dyy(V[z]) + Dyz(V[y]))
        +
        (gama*nu*nu)*(
            polmaginv*polmaginv*(-2*px*pz*(dxpx*px + dxpy*py + dxpz*pz)*(px*px*zdxvx + px*py*(zdxvy + zdyvx) + py*py*zdyvy + px*pz*(zdxvz + zdzvx) +
                py*pz*(zdyvz + zdzvy) + pz*pz*zdzvz) -
             2*py*pz*(dypx*px + dypy*py + dypz*pz)*(px*px*zdxvx + px*py*(zdxvy + zdyvx) + py*py*zdyvy + px*pz*(zdxvz + zdzvx) +
                py*pz*(zdyvz + zdzvy) + pz*pz*zdzvz) +
             1/3.0*(2*(dzpx*px + dzpy*py + dzpz*pz)*(px*px + py*py - 2*pz*pz)*
                (px*px*zdxvx + px*py*(zdxvy + zdyvx) + py*py*zdyvy + px*pz*(zdxvz + zdzvx) + py*pz*(zdyvz + zdzvy) +
                  pz*pz*zdzvz))) +
          polmaginv*((dxpz*px - (2*dzpx*px)/3. + dypz*py - (2*dzpy*py)/3. + dxpx*pz + dypy*pz + (4*dzpz*pz)/3.)*
              (px*px*zdxvx + px*py*(zdxvy + zdyvx) + py*py*zdyvy + px*pz*(zdxvz + zdzvx) + py*pz*(zdyvz + zdzvy) +
                pz*pz*zdzvz) + px*pz*(2*dxpx*px*zdxvx + px*px*zdxxvx + py*py*zdxxvy + pz*pz*zdxxvz + dxpy*px*(zdxvy + zdyvx) +
                dxpx*py*(zdxvy + zdyvx) + 2*dxpy*py*zdyvy + px*py*(zdxxvy + zdyxvx) + dxpz*px*(zdxvz + zdzvx) +
                dxpx*pz*(zdxvz + zdzvx) + dxpz*py*(zdyvz + zdzvy) + dxpy*pz*(zdyvz + zdzvy) + 2*dxpz*pz*zdzvz +
                px*pz*(zdxxvz + zdzxvx) + py*pz*(zdyxvz + zdzxvy)) +
             py*pz*(2*dypx*px*zdxvx + px*px*zdxyvx + dypy*px*(zdxvy + zdyvx) + dypx*py*(zdxvy + zdyvx) + 2*dypy*py*zdyvy +
                px*py*(zdxyvy + zdyyvx) + py*py*zdyyvy + pz*pz*zdyzvz + dypz*px*(zdxvz + zdzvx) + dypx*pz*(zdxvz + zdzvx) +
                dypz*py*(zdyvz + zdzvy) + dypy*pz*(zdyvz + zdzvy) + 2*dypz*pz*zdzvz + px*pz*(zdxyvz + zdzyvx) +
                py*pz*(zdyyvz + zdzyvy)) +
             1/3.0*((-px*px - py*py + 2*pz*pz)*
                (2*dzpx*px*zdxvx + px*px*zdxzvx + dzpy*px*(zdxvy + zdyvx) + dzpx*py*(zdxvy + zdyvx) + 2*dzpy*py*zdyvy +
                  px*py*(zdxzvy + zdyzvx) + py*py*zdyzvy + dzpz*px*(zdxvz + zdzvx) + dzpx*pz*(zdxvz + zdzvx) +
                  dzpz*py*(zdyvz + zdzvy) + dzpy*pz*(zdyvz + zdzvy) + 2*dzpz*pz*zdzvz + px*pz*(zdxzvz + zdzzvx) +
                  py*pz*(zdyzvz + zdzzvy) + pz*pz*zdzzvz)))
             ) 
        ;


        auto eq_sigma_xy=eta * (Dx(V[y]) + Dy(V[x])) +
                gama*nu*nu*polmaginv*(px*px*px*py*Dx(V[x]) + px*px*py*py*Dx(V[y]) + px*px*py*py*Dy(V[x]) + px*py*py*py*Dy(V[y]) +
                px*px*py*pz*Dx(V[z]) +  px*px*py*pz*Dz(V[x]) + px*py*py*pz*Dy(V[z]) + px*py*py*pz*Dz(V[y])+ px*py*pz*pz*Dz(V[z]));
        //xy rhs
        g_up[x]=-(zeta*delmu*px*py +  0.5*Hp3*px*px + 0.5*Hp3*py*py - 0.5*Hp1*px*pz - 0.5*Hp2*py*pz + (0.5*Hp3*px*px - 0.5*Hp3*py*py - 0.5*Hp1*px*pz + 0.5*Hp2*py*pz -
                        px*py*gama*lambda*delmu)*nu) - sigma[x][y];
        Particles.ghost_get<ExtForce>(SKIP_LABELLING);

        if(v_cl.rank()==0){
        std::cout << "Init of Velocity took " << tt.getwct() << " seconds." << std::endl;
        std::cout << "Calculate velocity (step t=" << t << ")" << std::endl;
        }
        tt.start();
        double V_err = 1, V_err_old,sum=0,sum1=0;
        int n = 0;
        int nmax = 30;
        int errctr, Vreset = 0;
        V_err = 1;
        n = 0;
        errctr = 0;
        eq_id vx, vy, vz;
        vx.setId(0);
        vy.setId(1);
        vz.setId(2);
        P = 0;
        V_t = 0;
        DCPSE_scheme<equations3d3Pxz, vector_type> Solver(Particles);
        Solver.impose(Stokes1, bulk, RHS[0], vx);
        Solver.impose(Stokes2, bulk, RHS[1], vy);
        Solver.impose(Stokes3, bulk, RHS[2], vz);
        Solver.impose(eq_sigma_xy, Boundary_up, g[x], vx);
        Solver.impose(V[y], Boundary_up, 0, vy);
        Solver.impose(V[z], Boundary_up, 0, vz);
        Solver.impose(V[x], Boundary_down, 0, vx);
        Solver.impose(V[y], Boundary_down, 0, vy);
        Solver.impose(V[z], Boundary_down, 0, vz);

        while (V_err >= V_err_eps && n <= nmax) {
            Particles.ghost_get<Pressure>(SKIP_LABELLING);
            RHS_bulk[x] = dV[x]+Dx(P);
            RHS_bulk[y] = dV[y]+Dy(P);
            RHS_bulk[z] = dV[z]+Dz(P);
            Particles.ghost_get<VRHS>(SKIP_LABELLING);
            Solver.reset_b();
            Solver.reset_x_ig();
            Solver.impose_b(bulk, RHS[0], vx);
            Solver.impose_b(bulk, RHS[1], vy);
            Solver.impose_b(bulk, RHS[2], vz);
            Solver.impose_b(Boundary_up, g[x], vx);
            Solver.impose_b(Boundary_up, 0, vy);
            Solver.impose_b(Boundary_up, 0, vz);
            Solver.impose_b(Boundary_down, 0, vx);
            Solver.impose_b(Boundary_down, 0, vy);
            Solver.impose_b(Boundary_down, 0, vz);
            Solver.impose_x_ig(bulk, V[x], vx);
            Solver.impose_x_ig(bulk, V[y], vy);
            Solver.impose_x_ig(bulk, V[z], vz);
            Solver.impose_x_ig(Boundary_up, V[x], vx);
            Solver.impose_x_ig(Boundary_up, V[y], vy);
            Solver.impose_x_ig(Boundary_up, V[z], vz);
            Solver.impose_x_ig(Boundary_down, V[x], vx);
            Solver.impose_x_ig(Boundary_down, V[y], vy);
            Solver.impose_x_ig(Boundary_down, V[z], vz);
            if(n==0){
                Solver.solve_with_solver(solverPetsc, V[x], V[y], V[z]);
            }
            else{
                Solver.solve_with_solver_ig_successive(solverPetsc, V[x], V[y], V[z]);
            }
            V_down[y]=0.0;
            Particles.ghost_get<Velocity>(SKIP_LABELLING);
            div = -(Dx(V[x]) + Dy(V[y])+Dz(V[z]));
            P = P + div;
            sum = 0;
            sum1 = 0;
            for (int j = 0; j < bulk.size(); j++) {
                auto p = bulk.get<0>(j);
                sum += (Particles.getProp<VT>(p)[0] - Particles.getProp<Velocity>(p)[0]) *
                       (Particles.getProp<VT>(p)[0] - Particles.getProp<Velocity>(p)[0]) +
                       (Particles.getProp<VT>(p)[1] - Particles.getProp<Velocity>(p)[1]) *
                       (Particles.getProp<VT>(p)[1] - Particles.getProp<Velocity>(p)[1]) +
                       (Particles.getProp<VT>(p)[2] - Particles.getProp<Velocity>(p)[2]) *
                       (Particles.getProp<VT>(p)[2] - Particles.getProp<Velocity>(p)[2]);
                sum1 += Particles.getProp<Velocity>(p)[0] * Particles.getProp<Velocity>(p)[0] +
                        Particles.getProp<Velocity>(p)[1] * Particles.getProp<Velocity>(p)[1]+
                        Particles.getProp<Velocity>(p)[2] * Particles.getProp<Velocity>(p)[2];
            }
            v_cl.sum(sum);
            v_cl.sum(sum1);
            v_cl.execute();
            sum = sqrt(sum);
            sum1 = sqrt(sum1);
            V_err_old = V_err;
            V_err = sum/sum1;
            if (V_err > V_err_old || abs(V_err_old - V_err) < 1e-8) {
                errctr++;
            } else {
                errctr = 0;
            }
            if (n > 5) {
                if (errctr > 3) {
                    if(v_cl.rank()==0){
                        std::cout << "CONVERGENCE LOOP BROKEN DUE TO INCREASE/VERY SLOW DECREASE IN ERROR : " <<V_err<< std::endl;
                    }
                    Vreset = 1;
                    break;
                } else {
                    Vreset = 0;
                }
            }
            V_t = V;
            n++;
           if (v_cl.rank() == 0) {
                std::cout << "Rel l2 cgs err in V: " << V_err << std::endl;
            }
        }
        tt.stop();
        double MaxRateOfChange=0;
        for (int j = 0; j < bulk.size(); j++) {
            auto p = bulk.get<0>(j);
            for (int i=0;i<3;i++){
                if(fabs((Particles.getProp<Polarization>(p)[i]-Particles.getProp<POLD>(p)[i])/(t-t_old))>MaxRateOfChange)
                {
                    MaxRateOfChange=fabs((Particles.getProp<Polarization>(p)[i]-Particles.getProp<POLD>(p)[i])/(t-t_old));
                }
            }
        }

        v_cl.max(MaxRateOfChange);
        v_cl.execute();
        if(v_cl.rank()==0)
        {
            std::cout<<"MaxRateOfChange: "<<MaxRateOfChange<<std::endl;
        }
        Pol_old=Pol;
        Particles.ghost_get<Velocity>(SKIP_LABELLING);
        u[x][x] = Dx(V[x]);
        u[x][y] = 0.5 * (Dx(V[y]) + Dy(V[x]));
        u[x][z] = 0.5 * (Dx(V[z]) + Dz(V[x]));
        u[y][x] = 0.5 * (Dy(V[x]) + Dx(V[y]));
        u[y][y] = Dy(V[y]);
        u[y][z] = 0.5 * (Dy(V[z]) + Dz(V[y]));
        u[z][x] = 0.5 * (Dz(V[x]) + Dx(V[z]));
        u[z][y] = 0.5 * (Dz(V[y]) + Dy(V[z]));
        u[z][z] = Dz(V[z]);

        if (v_cl.rank() == 0) {
            std::cout << "Rel l2 cgs err in V = " << V_err << " and took " << tt.getwct() << " seconds with " << n
                      << " iterations. dt is set to " << t-t_old
                      << std::endl;
        }

        t_old=t;
        W[x][x] = 0;
        W[x][y] = 0.5 * (Dy(V[x]) - Dx(V[y]));
        W[x][z] = 0.5 * (Dz(V[x]) - Dx(V[z]));
        W[y][x] = 0.5 * (Dx(V[y]) - Dy(V[x]));
        W[y][y] = 0;
        W[y][z] = 0.5 * (Dz(V[y]) - Dy(V[z]));
        W[z][x] = 0.5 * (Dx(V[z]) - Dz(V[x]));
        W[z][y] = 0.5 * (Dy(V[z]) - Dz(V[y]));
        W[z][z] = 0;


        if(ctr%wr_at==0 || ctr==wr_f){
        Particles.deleteGhost();
        Particles.write_frame("Polar3d", ctr,BINARY);
        Particles.ghost_get<0>();
        }
    }
};



int main(int argc, char* argv[])
{

    {   openfpm_init(&argc,&argv);
        timer tt2;
        tt2.start();
        //size_t grd_sz = int(std::atof(argv[1]));
        size_t grd_sz_z = 5;//int(std::atof(argv[2]));
        double tf = std::atof(argv[1]);
        dt = tf/std::atof(argv[2]);
        wr_f=int(std::atof(argv[2]));
        wr_at=int(std::atof(argv[3]));
        V_err_eps = double(std::atof(argv[4]));
        //size_t grd_sz=int(sqrt(8.0*std::atof(argv[2])/125.0)+1);
        size_t grd_sz=int(sqrt(0.64*std::atof(argv[2]))+1);
        
        double boxsize = 10.0 ;
        const size_t sz[3] = {grd_sz, grd_sz+1, grd_sz_z};
        std::cout<<"GridSize:"<<sz[0]<<"x"<<sz[1]<<"x"<<sz[2]<<std::endl;
        double spacing_p = 10.0 / (sz[0]);
        double spacing_np = 10.0 / (sz[1] - 1);
        Box<3, double> box({0, -4.0*spacing_np, 0}, {boxsize, boxsize+4.0*spacing_np, (grd_sz_z)*spacing_p});
        size_t bc[3] = {PERIODIC, NON_PERIODIC, PERIODIC};
        double Lx = 10;
        double Ly = 10;
        double Lz = (grd_sz_z)*spacing_p;
        double rCut = 3.9 * spacing_p;
        double rCut2 = 3.9 * spacing_p;
        int ord = 2;
        int ord2 = 2;
        double sampling_factor = 3.1;
        double sampling_factor2 = 1.9;
        Ghost<3, double> ghost(rCut2+spacing_np/8.0);
        Ghost<3,long int> gg(4);
        auto &v_cl = create_vcluster();


        vector_dist_ws<3, double, ActiveGel3d> Particles(0,box,bc,ghost);
        vector_dist<3, double, aggregate<double[3],double>> vd(0,box,bc,ghost);
        //creating a grid
        size_t gdsz[3]={sz[0],sz[1]+8,sz[2]};
        grid_type gd(Particles.getDecomposition(),gdsz,gg);

        Particles.setPropNames(PropNAMES);
        double x0, y0, z0, x1, y1, z1;
        x0 = 0;
        y0 = 0;
        z0 = 0;
        x1 = boxsize;
        y1 = boxsize;
        z1 = (grd_sz_z)*spacing_p;
        if (v_cl.rank() == 0) {
        std::cout<<"GridSize:"<<sz[0]<<"x"<<sz[1]<<"x"<<sz[2]<<std::endl;
        std::cout<<"Spacing Periodic:"<<spacing_p<<", Spacing Np:"<<spacing_np<<std::endl;}

        // create particles at grid nodes
        auto it = gd.getDomainIterator();

        while (it.isNext()) {
            auto lkey = it.get();
            auto key = it.getGKey(lkey);
            double y = key.get(1) * spacing_np;
            if(y<=y1 && y>=0.0){
                Particles.add();
                //Particles.getLastProp<GRIDID>() = lkey;
                double x = key.get(0) * spacing_p;
                Particles.getLastPos()[0] = x;
                Particles.getLastPos()[1] = y;

                double z = key.get(2) * spacing_p;
                Particles.getLastPos()[2] = z;
                if (y>10.0-spacing_np/8.0)
                {
                    Particles.getLastSubset(1);
                }
                else if (y<spacing_np/8.0)
                {
                    Particles.getLastSubset(2);
                }
                else
                {
                    Particles.getLastSubset(0);
                }

                // mirror particles
                if (y>10.0-(3*spacing_np+spacing_np/8.0) && y<10.0-spacing_np/8.0)
                {
                    //std::cout<<"domain is:"<<(p.getKey()>Particles.size_local())<<" and position is:("<<xp<<","<<yp<<","<<zp<<")"<<std::endl;
                    vd.add();
                    vd.getLastPos()[0]=x;
                    vd.getLastPos()[1]=(y)+2*(10.0-y);
                    vd.getLastPos()[2]=z;
                }
                else if (y<(3*spacing_np+spacing_np/8.0)&& y>0.0+spacing_np/8.0)
                {
                    vd.add();
                    vd.getLastPos()[0]=x;
                    vd.getLastPos()[1]=-y;
                    vd.getLastPos()[2]=z;
                }

                Particles.getLastProp<PARTICLENUMBER>() = key.get(0) + (key.get(1)*sz[0]) + (key.get(2)*sz[0]*sz[1]);   // individual number: usefull to keep track of particle in saved data
            }


            ++it;
        }

        Particles.map();
        vd.ghost_get<0,1>();
        Particles.ghost_get<Polarization,DPOL,DELMU,ExtForce>();

        auto Pol = getV<Polarization>(Particles);
        auto dPol= getV<DPOL>(Particles);
        auto Pol_old = getV<POLD>(Particles);
        auto g = getV<ExtForce>(Particles);
        auto V = getV<Velocity>(Particles);
        auto P = getV<Pressure>(Particles);
        auto delmu=getV<DELMU>(Particles);
        auto PolMag=getV<POLMAG>(Particles);
        g = 0;
        delmu = dmu;
        P = 0;
        V = 0;
        Pol_old=0;
        Particles.ghost_get<ExtForce>(SKIP_LABELLING);

        // set initial polarization
        auto it2 = Particles.getDomainIterator();

            while (it2.isNext()) {
                auto p = it2.get();
                Point<3, double> xp = Particles.getPos(p);
                Particles.getProp<0>(p)[x] = cos(M_PI * xp[y] / (2.0 * Ly));
                Particles.getProp<0>(p)[y] = sin(M_PI * xp[y] / (2.0 * Ly));
                Particles.getProp<0>(p)[z] = 0;
                if (xp[1]==y1)
                {
                    Particles.getProp<0>(p)[x] = 0;
                    Particles.getProp<0>(p)[y] = 1.0;
                    Particles.getProp<0>(p)[z] = 0;
                }
                ++it2;
            }

        vector_dist_subset<3, double, ActiveGel3d> Particles_bulk(Particles,0);
        vector_dist_subset<3, double, ActiveGel3d> Particles_boundary_up(Particles,1);
        vector_dist_subset<3, double, ActiveGel3d> Particles_boundary_dw(Particles,2);
        auto Pol_bulk = getV<Polarization>(Particles_bulk);
        auto & bulk = Particles_bulk.getIds();
        auto & Boundary_up = Particles_boundary_up.getIds();
        auto & Boundary_down = Particles_boundary_dw.getIds();

        // initialize DC PSE operators
        Derivative_x Dx(Particles, ord, rCut, sampling_factor, support_options::RADIUS);
        Derivative_y Dy(Particles, ord, rCut, sampling_factor, support_options::RADIUS);
        Derivative_z Dz(Particles, ord, rCut, sampling_factor, support_options::RADIUS);
        Derivative_xy Dxy(Particles, ord, rCut2, sampling_factor2, support_options::RADIUS);
        Derivative_yz Dyz(Particles, ord, rCut2, sampling_factor2, support_options::RADIUS);
        Derivative_xz Dxz(Particles, ord, rCut2, sampling_factor2, support_options::RADIUS);
        auto Dyx = Dxy;
        auto Dzy = Dyz;
        auto Dzx = Dxz;

        Derivative_xx Dxx(Particles, ord, rCut2, sampling_factor2,
                          support_options::RADIUS);//, Dxx2(Particles, ord2, rCut2, sampling_factor2, support_options::RADIUS);
        Derivative_yy Dyy(Particles, ord, rCut2, sampling_factor2,
                          support_options::RADIUS);
        Derivative_zz Dzz(Particles, ord, rCut2, sampling_factor2,
                          support_options::RADIUS);
        typedef boost::numeric::odeint::controlled_runge_kutta<boost::numeric::odeint::runge_kutta_cash_karp54< state_type_3d_ofp,double,state_type_3d_ofp,double,boost::numeric::odeint::vector_space_algebra_ofp>> controlled_stepper_type;
        boost::numeric::odeint::runge_kutta4< state_type_3d_ofp,double,state_type_3d_ofp,double,boost::numeric::odeint::vector_space_algebra_ofp> rk4;
        vectorGlobal=(void *) &Particles;
        vectorGlobal_bulk=(void *) &Particles_bulk;
        vectorGlobal_boundary_up=(void *) &Particles_boundary_up;
        vectorGlobal_boundary_down=(void *) &Particles_boundary_dw;
        global_gd=(void *) &gd;
        global_vd=(void *) &vd;
        inte_type2m inte2m(vd,gd);
        inte_type2p inte2p(Particles,gd);



        PolarEv<Derivative_x,Derivative_y,Derivative_z,Derivative_xx,Derivative_xy,Derivative_xz,Derivative_yy,Derivative_yz,Derivative_zz> System(Dx,Dy,Dz,Dxx,Dxy,Dxz,Dyy,Dyz,Dzz);
        CalcVelocity<Derivative_x,Derivative_y,Derivative_z,Derivative_xx,Derivative_xy,Derivative_xz,Derivative_yy,Derivative_yz,Derivative_zz,inte_type2m,inte_type2p> CalcVelocityObserver(Dx,Dy,Dz,Dxx,Dxy,Dxz,Dyy,Dyz,Dzz,inte2m,inte2p);


        state_type_3d_ofp tPol;

        tPol.data.get<0>()=Pol[x];
        tPol.data.get<1>()=Pol[y];
        tPol.data.get<2>()=Pol[z];
        dPol=Pol;
    
        PolMag=Pol[x]*Pol[x]+Pol[y]*Pol[y]+Pol[z]*Pol[z];
        Particles.ghost_get<Polarization>(SKIP_LABELLING);
        double tim = 0; // intermediate time steps
        std::vector<double> inter_times;
        std::ofstream myFile("PolMag");
        myFile <<argv[1]<<","<<argv[2]<<","<<argv[3]<<","<<argv[4]<<","<< "\n";
        myFile.close();
        size_t steps=boost::numeric::odeint::integrate_const(rk4,System,tPol,0.0,tf,dt,CalcVelocityObserver);


        std::cout << "Time steps: " << steps << std::endl;

        Pol_bulk[x]=tPol.data.get<0>();
        Pol_bulk[y]=tPol.data.get<1>();
        Pol_bulk[z]=tPol.data.get<2>();


        Particles.deleteGhost();
        Particles.write_frame("anaHyb",wr_f,BINARY);
        Particles.save("HybEND");
        Dx.deallocate(Particles);
        Dy.deallocate(Particles);
        Dz.deallocate(Particles);
        Dxy.deallocate(Particles);
        Dxz.deallocate(Particles);
        Dyz.deallocate(Particles);
        Dxx.deallocate(Particles);
        Dyy.deallocate(Particles);
        Dzz.deallocate(Particles);
        Particles.deleteGhost();
        tt2.stop();
        if (v_cl.rank() == 0) {
            std::cout << "The simulation took " << tt2.getcputime() << "(CPU) ------ " << tt2.getwct()
                      << "(Wall) Seconds.";
        }
    }
    openfpm_finalize();

}