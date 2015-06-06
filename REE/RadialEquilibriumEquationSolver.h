#include <vector> 

struct data {
	std::vector<double> r;
	std::vector<double> z;
};

class REESolver
{
public:
	REESolver();
	~REESolver();
	bool runREESolver();
	bool setData();
private:
	int _nStation;
	int _nStream;
	int _indexLE;
	int _indexTE;
	double _gamma;
	double _gasConst;
	double _cp;
	double _tolPsi;
	double _tolRho;
	double _tolRHS;
	double _convCritRho;
	double _convCritRHS;
	double _convCritPsi;
	int _maxIt;
	double _rho0In;
	double _p0In;
	double _T0In;
	double _rHub;
	double _rShd;
	double _chord;
	double _rCtIn;
	double _rCtOut;
	double _mfr;
	double _omega;
	double _lossesTE;
	double _rothalpy;
	std::vector<double> _losses;
	std::vector<double> _beta;
	data _mesh;
	std::vector<std::vector<double> > _psi;
	std::vector<std::vector<double> > _Cz;
	std::vector<std::vector<double> > _Cr;
	std::vector<std::vector<double> > _rCt;
	std::vector<std::vector<double> > _rho;
	std::vector<std::vector<double> > _h0;
	std::vector<std::vector<double> > _p0;
	std::vector<std::vector<double> > _entropy;
	std::vector<std::vector<double> > _temp;
	std::vector<std::vector<double> > _RHS;

	bool initialize();
	bool initializeField();
	bool computeVelocityField();
	bool computeThermoVariables();
	bool computeRHS();
	bool updateStreamField();
	double computeMaxChangeInRHS();
	double computeMaxChangeInDensity();
	bool writeReport();
};