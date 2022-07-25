/*
    Copyright (C) 2022 Giacomo Rosilho de Souza

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include "MonteCarlo.h"
#include "ProblemsList.h"
#include "SSA.h"
#include <sys/stat.h>

static inline void loadBar(int x, int n, double start, double end);

MonteCarlo::MonteCarlo(Parameters param_)
:param(param_)
{
    MCiter = param.get_MCiter();
    X.assign(MCiter, VectorXd::Zero(param.get_num_species()));
    mean = ArrayXd::Zero(param.get_num_species());
    std_dev = ArrayXd::Zero(param.get_num_species());
}

MonteCarlo::~MonteCarlo()
{
}

void MonteCarlo::run()
{    
    #ifdef _OPENMP
    unsigned int num_threads = omp_get_max_threads();
    #else
    unsigned int num_threads = 1;
    #endif

    elapsed_time = get_cpu_time()/num_threads;
    
    vector<bool> succeed(MCiter,false);
    vector<Real> tau_mean_vec(MCiter);
    vector<Real> s_mean_vec(MCiter);
    vector<Real> damping_mean_vec(MCiter);
    vector<Real> rho_Nit_mean_vec(MCiter);
    
    #pragma omp parallel 
    {
        ChemicalSystem* cs;
        Solver* sol;
        param.init_Solver_ChemicalSystem(sol,cs);
        
        #pragma omp for schedule(dynamic, 10)
        for(unsigned int i=0; i<MCiter; i++)
        {
            if(sol->run())
            {
                X[i] = sol->get_sol();
                tau_mean_vec[i] = sol->get_tau_mean();
                s_mean_vec[i] = sol->get_s_mean();
                rho_Nit_mean_vec[i] = sol->get_rho_or_Nit_mean();
                damping_mean_vec[i] = sol->get_damping_mean();
                succeed[i] = true;
            }
            loadBar(i, MCiter, elapsed_time, get_cpu_time()/num_threads);
        }
        
        delete sol;
        delete cs;
    }
    elapsed_time = get_cpu_time()/num_threads-elapsed_time;
    
    process_data(succeed,tau_mean_vec,s_mean_vec,rho_Nit_mean_vec,damping_mean_vec);
    write_data(succeed);
    print_data(param.get_solver(), param.get_post_proc(), param.get_Newton_tol());
    
    if(param.get_refsol()!=string(""))
        compute_print_errors();
}

void MonteCarlo::process_data(vector<bool>& succeed, vector<Real>& tau_mean_vec,
                 vector<Real>& s_mean_vec, vector<Real>& rho_Nit_mean_vec, 
                 vector<Real>& damping_mean_vec)
{    
    cout<<"Processing data..."<<endl<<flush;
    
    tau_mean=0;
    s_mean=0;
    rho_Nit_mean=0.;
    damping_mean=0.;
    succeed_MCiter = 0;
    for(unsigned int i=0;i<MCiter;i++)
        if(succeed[i])
        {
            mean += X[i].array();
            tau_mean += tau_mean_vec[i];
            s_mean += s_mean_vec[i];
            rho_Nit_mean += rho_Nit_mean_vec[i];
            damping_mean += damping_mean_vec[i];
            succeed_MCiter++;
        }
    mean /= succeed_MCiter;
    tau_mean /= succeed_MCiter;
    s_mean /= succeed_MCiter;
    rho_Nit_mean /= succeed_MCiter;
    damping_mean /= succeed_MCiter;
    
    for(unsigned int i=0;i<MCiter;i++)
        if(succeed[i])
            std_dev += (X[i].array()-mean)*(X[i].array()-mean); // operations are element-wise, are we are using array instead of vector
    std_dev /= (succeed_MCiter-1.); 
    std_dev = std_dev.array().sqrt();
       
    printf("\033[F\033[J");
}

void MonteCarlo::print_data(string solver_name, bool post_proc, Real Ntol)
{
    if(param.get_refsol()!=string(""))
        cout<<endl<<"----------------- SOLUTION DATA -----------------"<<endl<<endl;
    
    print_solver_info(elapsed_time,MCiter,succeed_MCiter,solver_name,
                      tau_mean,s_mean,damping_mean,rho_Nit_mean,post_proc,Ntol);
     
    print_statistics(mean, std_dev);
}

void MonteCarlo::print_solver_info(Real time, unsigned int iter, unsigned int succ_iter, 
                 string solver_name, Real tau_mean, Real s_mean, Real damping_mean, 
                 Real rho_Nit_mean, bool postproc, Real Ntol)
{
    cout.unsetf(ios::fixed | ios::scientific);
    cout<<setfill(' ');
    cout<<setprecision(4)<<scientific;
    
    cout<<"-------------- Solver info -------------"<<endl;
    cout<<setw(20)<<left<<"Solver: "<<solver_name<<endl;
    cout<<setw(20)<<left<<"Mean nsteps : "<<param.get_final_time()/tau_mean<<endl;
    cout<<setw(20)<<left<<"Mean tau : "<<tau_mean<<endl;
    cout<<setw(20)<<left<<"Mean s : "<<s_mean<<endl;
    cout<<setw(21)<<left<<"Mean \u03B5 : "<<damping_mean<<endl;
    cout<<setw(21)<<left<<"Mean \u03C1/Nit : "<<rho_Nit_mean<<endl;
    cout<<setw(20)<<left<<"Newton tol : "<<Ntol<<endl;
    cout<<setw(20)<<left<<"Postprocessing: "<<(postproc? "yes":"no")<<endl;
    cout<<setw(20)<<left<<"Time to solution: "<<time<<" sec"<<endl;
    cout<<setw(20)<<left<<"Monte Carlo iter: "<<iter<<endl;
    cout<<setw(20)<<left<<"Successful iter: "<<succ_iter<<endl;
    
    cout.unsetf(ios::fixed | ios::scientific);
    cout<<setw(20)<<left<<"Succeed ratio: "<<100.*succ_iter/iter<<" %"<<endl;
    #ifdef _OPENMP
    unsigned int num_threads = omp_get_max_threads();
    #else
    unsigned int num_threads = 1;
    #endif
    cout<<setw(20)<<left<<"Number of threads: "<<num_threads<<endl;
    cout<<"----------------------------------------"<<endl;
}

void MonteCarlo::print_statistics(const ArrayXd& avg, const ArrayXd& std)
{
    unsigned int prec = 4;
    unsigned int num_size = prec+6;
    unsigned int cell_size = num_size+1;
    unsigned int N = avg.size();
    unsigned int table_width = 6+(cell_size+1)*N;
    
    unsigned int short_width = (table_width-12);
    
    cout<<setprecision(prec)<<scientific;
    
    
    cout<<setfill('-')<<setw(short_width/2)<<"-"<<" Statistics "
        <<setw(short_width/2+short_width%2)<<"-"<<endl;

    cout<<"|Avg |";
    for(unsigned int i=0;i<N;i++)
        cout<<setfill(' ')<<setw(cell_size)<<right<<avg(i)<<"|";
    cout<<endl;
    
    cout<<setfill('-')<<setw(table_width)<<"-"<<endl;
    cout<<"|Std |";
    for(unsigned int i=0;i<N;i++)
        cout<<setfill(' ')<<setw(cell_size)<<right<<std(i)<<"|";
    cout<<endl;
    cout<<setfill('-')<<setw(table_width)<<"-"<<endl;
}

void MonteCarlo::write_data(vector<bool>& succeed)
{
    cout<<"Writing data..."<<endl<<flush;
    
    // write all samples in a binary file
    unsigned int N = X[0].size();
    
    string filename = param.get_filename()+string("_data.bin");
    remove(filename.c_str());
    ofstream ofile(filename, ios::out | ios::binary);
    ofile.seekp(0);

    for(unsigned int i=0;i<MCiter;i++)
        if(succeed[i])
            ofile.write((char*)&X[i](0), N*sizeof(double));
    
    ofile.close();
    
    // write statistics in a csv file
    filename = param.get_filename()+string("_statistics.csv");
    remove(filename.c_str());
    ofile.open(filename, ios::out);
    
    ofile<<setprecision(16);
    ofile<<"mean, std";
    for(unsigned int i=0;i<mean.size();i++)
        ofile<<endl<<mean(i)<<", "<<std_dev(i);
    
    ofile.close();
    
    // same as before but in a binary file
    filename = param.get_filename()+string("_statistics.bin");
    remove(filename.c_str());
    ofile.open(filename, ios::out | ios::binary);
    ofile.seekp(0);
    
    for(unsigned int i=0;i<mean.size();i++)
    {
        ofile.write((char*)&mean(i),sizeof(double));
        ofile.write((char*)&std_dev(i),sizeof(double));
    }
    
    ofile.close();
    
    // now write solver info in a csv file
    filename = param.get_filename()+string("_MCinfo.csv");
    remove(filename.c_str());
    ofile.open(filename, ios::out);
    
    ofile<<setprecision(16);
    ofile<<"time, MCiter, succeed_MCiter, solver, tau_mean, s_mean, damping_mean, rho_Nit_mean, postproc, newton_tol"<<endl;
    ofile<<elapsed_time<<", "<<MCiter<<", "<<succeed_MCiter
         <<", "<<param.get_solver()<<", "<<tau_mean<<", "<<s_mean<<", "<<damping_mean<<", "<<rho_Nit_mean<<", "<<param.get_post_proc()<<", "<<param.get_Newton_tol();
    
    ofile.close();
    
    // same as before but in binary file
    filename = param.get_filename()+string("_MCinfo.bin");
    remove(filename.c_str());
    ofile.open(filename, ios::out | ios::binary);
    ofile.seekp(0);
    
    unsigned int solv_string_size =param.get_solver().length()+1;
    char solver_name[solv_string_size];
    strcpy(solver_name,param.get_solver().c_str());
    bool postproc = param.get_post_proc();
    Real Ntol = param.get_Newton_tol();
    
    ofile.write((char*)&elapsed_time,sizeof(double));
    ofile.write((char*)&MCiter,sizeof(unsigned int));
    ofile.write((char*)&succeed_MCiter,sizeof(unsigned int));
    ofile.write((char*)&solv_string_size,sizeof(unsigned int)); //needed when reading the same file
    ofile.write((char*)&solver_name,sizeof(char)*solv_string_size);
    ofile.write((char*)&tau_mean,sizeof(double));
    ofile.write((char*)&s_mean,sizeof(double));
    ofile.write((char*)&damping_mean,sizeof(double));
    ofile.write((char*)&rho_Nit_mean,sizeof(double));
    ofile.write((char*)&postproc,sizeof(bool));
    ofile.write((char*)&Ntol,sizeof(double));
    
    ofile.close();
    
    printf("\033[F\033[J");
}

void MonteCarlo::compute_print_errors()
{
    cout<<"Importing reference solution data..."<<endl<<flush;
    
    //import reference solution
    vector<VectorXd> refX;
    import_data(param.get_refsol(),refX);
    
    printf("\033[F\033[J");
    cout<<"Reading reference solution solver info..."<<endl<<flush;
    
    // reading SOLVER INFO of ref solution
    Real elapsed_time_ref;
    unsigned int MCiter_ref, succeed_MCiter_ref;
    Real tau_mean,s_mean,damping_mean,rho_Nit_mean,Ntol;
    bool postproc_ref;
    string solver_name_ref;
    
    read_solver_info(param.get_refsol(),elapsed_time_ref,MCiter_ref,succeed_MCiter_ref,
            solver_name_ref,tau_mean,s_mean,damping_mean,rho_Nit_mean,postproc_ref,Ntol);
        
    printf("\033[F\033[J");
    cout<<"Reading reference solution statistics..."<<endl<<flush;
    
    // reading statistics of ref sol
    ArrayXd mean_ref, std_dev_ref;
    read_statistics(param.get_refsol(),mean_ref,std_dev_ref);
        
    printf("\033[F\033[J");
    cout<<"Computing density distance area..."<<endl<<flush;
    
    vector<unsigned int> n_bins;
    ArrayXd dda;
    compute_density_distance_area(n_bins,dda,X,refX);
        
    printf("\033[F\033[J");
    cout<<"Estimating self distances..."<<endl<<flush;
    
    ArrayXd sd1, sd2;
    estimate_self_distances(n_bins, X.size(), refX.size(), sd1, sd2);
    
    
    printf("\033[F\033[J");
    
    // Printing info of reference solution
    cout<<endl<<"------------ REFERENCE SOLUTION DATA ------------"<<endl<<endl;
    print_solver_info(elapsed_time_ref,MCiter_ref,succeed_MCiter_ref,
            solver_name_ref,tau_mean,s_mean,damping_mean,rho_Nit_mean,postproc_ref,Ntol);
    print_statistics(mean_ref, std_dev_ref);
     
    
    // printing errors
    cout<<endl<<"-------------------- ERRORS --------------------"<<endl<<endl;
    
    // print error of mean and std
    unsigned int N = param.get_num_species();
    ArrayXd mean_err = (mean-mean_ref)/mean_ref;
    ArrayXd std_dev_err = (std_dev-std_dev_ref)/std_dev_ref;
    mean_err = mean_err.abs();
    std_dev_err = std_dev_err.abs();
    unsigned int prec = 4;
    unsigned int num_size = prec+6;
    unsigned int cell_size = num_size+1;
    unsigned int table_width = 6+(cell_size+1)*N;
    unsigned int short_width = (table_width-12);
    
    cout<<setfill('-')<<setw(short_width/2-2)<<"-"<<" Relative Error "
        <<setw(short_width/2+short_width%2-2)<<"-"<<endl;

    cout<<"|Avg |";
    for(unsigned int i=0;i<N;i++)
        cout<<setfill(' ')<<setw(cell_size)<<right<<mean_err(i)<<"|";
    cout<<endl;
    
    cout<<setfill('-')<<setw(table_width)<<"-"<<endl;
    cout<<"|Std |";
    for(unsigned int i=0;i<N;i++)
        cout<<setfill(' ')<<setw(cell_size)<<right<<std_dev_err(i)<<"|";
    cout<<endl;
    cout<<setfill('-')<<setw(table_width)<<"-"<<endl;
    
    string filename = param.get_filename()+string("_errors.csv");
    remove(filename.c_str());
    ofstream ofile(filename, ios::out);
    for(unsigned int i=0;i<mean_err.size();i++)
    {
        ofile<<mean_err(i)<<", "<<std_dev_err(i)<<endl;
    }
    ofile.close();
    
    // printing the nbins, self distances, density distance area
    cout<<endl;
    cout<<setfill('-')<<setw(short_width/2-6)<<"-"<<" Densisty distance area "
        <<setw(short_width/2+short_width%2-6)<<"-"<<endl;

    cout<<"|nbin|";
    for(unsigned int i=0;i<N;i++)
        cout<<setfill(' ')<<setw(cell_size)<<right<<n_bins[i]<<"|";
    cout<<endl;
    cout<<setfill('-')<<setw(table_width)<<"-"<<endl;
    
    cout<<"|SDX |";
    for(unsigned int i=0;i<N;i++)
        cout<<setfill(' ')<<setw(cell_size)<<right<<sd1(i)<<"|";
    cout<<endl;
    cout<<setfill('-')<<setw(table_width)<<"-"<<endl;
    
    cout<<"|SDrX|";
    for(unsigned int i=0;i<N;i++)
        cout<<setfill(' ')<<setw(cell_size)<<right<<sd2(i)<<"|";
    cout<<endl;
    cout<<setfill('-')<<setw(table_width)<<"-"<<endl;
    
    cout<<"|DDA |";
    for(unsigned int i=0;i<N;i++)
        cout<<setfill(' ')<<setw(cell_size)<<right<<dda(i)<<"|";
    cout<<endl;
    cout<<setfill('-')<<setw(table_width)<<"-"<<endl;
}

void MonteCarlo::compute_density_distance_area(vector<unsigned int>& n_bins, ArrayXd& dda,
        vector<VectorXd>& X1, vector<VectorXd>& X2)
{
    unsigned int N = param.get_num_species();
    n_bins.resize(N,param.get_n_bins()); // default n bins
    dda.resize(N);
    
    for(unsigned int j=0;j<N;j++)
    {
        Real min_val = X1[0](j);
        Real max_val = X1[0](j);
        
        #pragma omp parallel reduction(min:min_val) reduction(max:max_val)
        {
        #pragma omp for
        for(unsigned int i=0;i<X1.size();i++)
        {
            min_val = min(min_val,X1[i](j));
            max_val = max(max_val,X1[i](j));
        }
        #pragma omp for
        for(unsigned int i=0;i<X2.size();i++)
        {
            min_val = min(min_val,X2[i](j));
            max_val = max(max_val,X2[i](j));
        }
        }
        
        // adaptive bins selection
        min_val = floor(min_val);
        if(max_val != ceil(max_val))
            max_val = ceil(max_val);
        else
            max_val++;
        if((int)round(max_val-min_val)%2==1)
            max_val++;
        
        if(n_bins[j]==0)
        {
            n_bins[j] = max_val-min_val; //bins of size 1 integer
            //but many variables increase or decrease their value by two, so its better to have bins of size 2
            n_bins[j] = n_bins[j]/2;
        }// else use the default (not zero) value
        
        dda(j)=0.;
        Real nX1, nX2;
        Real bin_min, bin_max;
        for(unsigned int b=0;b<n_bins[j];b++)
        {
            nX1=0.;
            nX2=0.;
            bin_min = min_val + b*((max_val-min_val)/n_bins[j]);
            bin_max = min_val + (b+1.)*((max_val-min_val)/n_bins[j]);
            
            #pragma omp parallel reduction(+:nX1,nX2)
            {
            #pragma omp for
            for(unsigned int i=0;i<X1.size();i++)
                if(bin_min<=X1[i](j) && X1[i](j)<bin_max)
                    nX1++;
            #pragma omp for
            for(unsigned int i=0;i<X2.size();i++)
                if(bin_min<=X2[i](j) && X2[i](j)<bin_max)
                    nX2++;
            }

            dda(j) += abs(nX1/X1.size()-nX2/X2.size());
        }
    }
    
}

void MonteCarlo::estimate_self_distances(vector<unsigned int>& n_bins, unsigned int N1, unsigned int N2,
        ArrayXd& sd1, ArrayXd& sd2)
{
    unsigned int N = n_bins.size();
    sd1.resize(N);
    sd2.resize(N);
    for(unsigned int i=0;i<N;i++)
    {
        sd1(i) = sqrt(4.*n_bins[i]/3.141592/N1);
        sd2(i) = sqrt(4.*n_bins[i]/3.141592/N2);
    }
}

void MonteCarlo::import_process_data()
{
    string solver_name;
    bool post_proc;
    Real Ntol;
    
    import_data(param.get_filename(),X);
    read_solver_info(param.get_filename(), elapsed_time, MCiter, succeed_MCiter, 
            solver_name,tau_mean,s_mean,damping_mean,rho_Nit_mean,post_proc,Ntol);
    read_statistics(param.get_filename(),mean,std_dev);
    
    print_data(solver_name, post_proc, Ntol);
    
    compute_print_errors();
}

void MonteCarlo::import_data(string filename, vector<VectorXd>& Y)
{
    //import solution data
    filename = filename + string("_data.bin");

    struct stat st;
    stat(filename.c_str(), &st);
    unsigned int file_size = st.st_size;
    
    unsigned int N = param.get_num_species();
    unsigned int samples = file_size/sizeof(double)/N;
    Y.resize(samples, VectorXd::Zero(N));
    
    ifstream ifile(filename, ios::in | ios::binary);
    for(unsigned int i=0;i<samples;i++)
            ifile.read((char*)&Y[i](0), N*sizeof(double));
    ifile.close();
    
}

void MonteCarlo::read_solver_info(string filename, Real& time, unsigned int& iter,
                 unsigned int& succeed_iter, string& solver_name_str,
                 Real& tau_mean, Real& s_mean, Real& damping_mean,Real& rho_Nit_mean, 
                 bool& post_proc, Real& Ntol)
{
    unsigned int solv_string_size;
    filename = filename+string("_MCinfo.bin");
    ifstream ifile(filename, ios::in | ios::binary);
    
    ifile.read((char*)&time,sizeof(double));
    ifile.read((char*)&iter,sizeof(unsigned int));
    ifile.read((char*)&succeed_iter,sizeof(unsigned int));
    ifile.read((char*)&solv_string_size,sizeof(unsigned int)); 
    char solver_name_char[solv_string_size];
    ifile.read((char*)&solver_name_char,sizeof(char)*solv_string_size);
    solver_name_str = solver_name_char;
    ifile.read((char*)&tau_mean,sizeof(double));
    ifile.read((char*)&s_mean,sizeof(double));
    ifile.read((char*)&damping_mean,sizeof(double));
    ifile.read((char*)&rho_Nit_mean,sizeof(double));
    ifile.read((char*)&post_proc,sizeof(bool));
    ifile.read((char*)&Ntol,sizeof(double));
    
    ifile.close();
}

void MonteCarlo::read_statistics(string filename, ArrayXd& avg, ArrayXd& std)
{
    unsigned int N = param.get_num_species();
    avg.resize(N);
    std.resize(N);
    filename = filename+string("_statistics.bin");
    ifstream ifile(filename, ios::in | ios::binary);
    
    for(unsigned int i=0;i<N;i++)
    {
        ifile.read((char*)&avg(i),sizeof(double));
        ifile.read((char*)&std(i),sizeof(double));
    }
    
    ifile.close();
}

Real MonteCarlo::get_cpu_time()
{
    return (Real)clock() / CLOCKS_PER_SEC;
}

// Process has done i out of n rounds,
// and we want a bar of width w and resolution r.
static inline void loadBar(int i, int n, double start, double end)
{
    // Only update r times.
    int r=n/100;
    if(r==0)
        r=1;
    if ( i % r != 0 ) 
        return;
 
    // Calculate the ratio of complete-to-incomplete.
    int w = 60;
    float ratio = i/(float)n;
    int   c     = ratio * w;
    float time = (end-start)/60;
    float tend = time/ratio-time;
 
    cout<<"--------------------------------------------------------------"<<endl<<flush;
    
    cout<<"Complete: "<<(int)(ratio*100)<<" %"<<endl;
    cout<<setprecision(1);
    cout<<"Elapsed time: "<<time<<" minutes."<<endl;
    cout<<"Time to end : "<<tend<<" minutes."<<endl;   
    
    // Show the load bar.
    cout<<"[";
    for (int x=0; x<c; x++)
       printf("=");
 
    for (int x=c; x<w; x++)
       printf(" ");
    cout<<"]"<<endl;
    cout<<"--------------------------------------------------------------"<<endl<<flush;
    
    printf("\033[F\033[J\033[F\033[J\033[F\033[J\033[F\033[J\033[F\033[J\033[F\033[J");
}

