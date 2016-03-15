#include <iostream>
#include <cstdlib>
#include <Vector.h>
#include <string.h>
#include <random>

std::mt19937 make_seeded_engine() {
    std::random_device r;
    std::seed_seq seed{r(), r(), r(), r(), r(), r(), r(), r()};
    return std::mt19937(seed);
}


class mRandom {
    public:
        mRandom(long int seed) {
            engine = make_seeded_engine() ; 
        }
        
        double getUniform() { 
            return dist_real (engine) ; 
        }
        double getNormal() {
            return dist_normal(engine) ; 
        }
        
        
        ~mRandom() = default ; 
    private:
        std::mt19937 engine ; 
        std::uniform_real_distribution<> dist_real {0, 1};
        std::normal_distribution<> dist_normal {0,1};

} ;


void help() {
    std::cout << "NP DIAM EPS RC\n" ;
}

mVector getRandomOrientation(mRandom * rg) {
     mVector v (0,0,0) ; 
    v(0)=rg->getNormal() ; 
    v(1)=rg->getNormal() ;
    v(2)=rg->getNormal() ;
    v=v/v.norm() ;
    return v ; 
}

double Energy (double EPS, double RM, double r) {
    return EPS*(std::pow(r/RM,-12)-2*(std::pow(r/RM,-6))) ; 
}


void DLAFactory(mRandom * mrg,int n, const double DIAM, const double EPS, const double RC) {
    
    mVector * AGG = new mVector [n] ; 
    AGG[0].zeros()  ; 
    const double STEPSIZE = DIAM/10 ; 
    const double MAXSTEP = 1000 ;
    const double BETA = 1 ; 
    double RM = 0 ; 
    int countDone = 1 ; 
    double gyrationRadius = DIAM/2 ; 
    do {
        double startingPosition = 2 * gyrationRadius ;
        mVector walker = startingPosition * getRandomOrientation(mrg) ; 
        bool contact = false ;
        bool lost = false ; 
        int test = 0 ;
        // std::cout << "Starting at "<< walker.norm() << "\n"; 
        do {
            mVector oldWalker = walker ; 
 
            do { 
                walker = oldWalker+ STEPSIZE* getRandomOrientation(mrg) ; 
            
                RM  = RC*gyrationRadius;
                double oldpsi=Energy(EPS,RM, oldWalker.norm()) ; 
                double newpsi=Energy(EPS,RM, walker.norm()) ; 
                double p = std::exp(-BETA*(newpsi-oldpsi)) ; 
                if (mrg->getUniform()<p) {
                    //std::cout << oldWalker.norm()<<"\t" << walker.norm() << "\t" << newpsi << "\t"<<  oldpsi << "\t" << p << "\n" ; 
                    break ;
                } 
            } while (true) ; 
                
          //  std::cout << walker.norm() << "\t" << RM<< "\n" ; 
            
            // Find Contact 
            for (int i = 0 ; i < countDone ; i++) {
                double d = (AGG[i]-walker).norm() ;
                if (d<DIAM) { 
                    contact=true ;
                    mVector pos=walker-AGG[i] ;
                    pos=pos*(DIAM/pos.norm()) ; 
                    AGG[countDone] = AGG[i]+pos ; 
                    countDone++;
                    break ;
                }
            }
            lost = (walker.norm()>2*startingPosition) ? true : false ; 
            test++ ; 
        } while ((!contact)&&(!lost)&&(test<MAXSTEP)) ; 
        //      find position of center of mass  :
        mVector g(0,0,0) ; 
        for (int i = 0 ; i < countDone ; i ++) g+=AGG[i] ;
        g=g/countDone ; 
        for (int i = 0 ; i < countDone ; i ++) AGG[i]-=g ;
        // Find Structure Size : 
        double sum =0;
        for (int i = 0 ; i < countDone ; i ++) sum+=AGG[i].squaredNorm() ;
        gyrationRadius=sqrt(sum/countDone)+DIAM*0.5 ;
        
    } while (countDone<n) ;
   for (int i = 0 ; i < n ; i++) std::cout << AGG[i] << "\n" ; 
   std::cout << gyrationRadius << "\n" ; 
    delete [] AGG ; 
}



int main ( int argc, char **argv ) {
    int n = 0 ;
    double diam = 0.0; 
    double eps = 0.0; 
    double rc = 0.0; 
    
    if (argc != 5) {
        help() ; 
        exit(1) ; 
    }
    std::stringstream ss;
    ss<<argv[1] ; 
    if (!(ss >> n))    std::cerr << "Invalid number " << argv[1] << '\n';
   ss.clear() ;  ss << argv[2];
    if (!(ss >> diam))    std::cerr << "Invalid number " << argv[2] << '\n';
   ss.clear() ; ss << argv[3];
    if (!(ss >> eps))    std::cerr << "Invalid number " << argv[3] << '\n';
    ss.clear() ; ss << argv[4];
    if (!(ss >> rc))    std::cerr << "Invalid number " << argv[4] << '\n';
    // NP DIAM EPS RC
    
    mRandom * mrg = new mRandom(0) ; 
    
    
    DLAFactory(mrg,n,diam,eps,rc) ; 
    
    
    delete mrg ; 
    exit(0) ; 
}
