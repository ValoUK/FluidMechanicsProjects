#include "RadialEquilibriumEquationSolver.h"
#include <fstream>
#include <iostream>
#include <algorithm>    // std::max, std::fill
#define _USE_MATH_DEFINES
#include <math.h> //log ...
#include <cmath> //M_PI

REESolver::REESolver()
{
	initialize();
}

REESolver::~REESolver()
{

}

bool REESolver::initialize()
{
	_gasConst = 287.1;
	_gamma = 1.4;

	_nStream = 11;
	_nStation = _nStream * 5 - 4;
	_indexLE = 20;
	_indexTE = 30;
	_maxIt = 200;
	_tolPsi = 1.e-5;
	_tolRho = 1.e-3;
	_tolRHS = 1.e-2;
	_convCritRho = 0.;
	_convCritRHS = 0.; 
	_convCritPsi = 0.;

	_rCtIn = 39.3;
	_rCtOut = 117.8;
	_chord = _rCtOut - _rCtIn;
	_rHub = 0.45;
	_rShd = 0.5;
	_T0In = 288;
	_rho0In = 1.5;
	_lossesTE = 0.03;
	_omega = 6000.; _omega = _omega*M_PI / 30.;
	_mfr = 30.4;
	_cp = _gamma*_gasConst / (_gamma - 1.);
	_p0In = _rho0In*_gasConst*_T0In;

	return true;
}

bool REESolver::setData() 
{
	
	return true; 
}

bool REESolver::runREESolver()
{
	int count(0);

	if(!setData()) {
		std::cout << "Error occur while setting data" << std::endl;
		return false;
	}
	initializeField();
	while ((_convCritRHS > _tolRHS) || (_convCritRho > _tolRho) ||
		count < 3) {
		computeVelocityField();
		computeThermoVariables();
		computeRHS();
		updateStreamField();
		count++;
		if (count>_maxIt) { 
			std::cout << "Maximum number of iterations reached in runREESolver()" << std::endl;
			return false;
		}
	}
	if(!writeReport())
		return false;
	return true;
}

bool REESolver::initializeField()
{
	_mesh.r.resize(_nStream, -1.);
	_mesh.z.resize(_nStation, -1.);
	_beta.resize(_nStream,-1.);
	std::vector<double> tmp(int(_mesh.r.size()), -1.);
	_psi.resize(int(_mesh.z.size()), tmp);
	_Cz.resize(int(_mesh.z.size()), tmp);
	_Cr.resize(int(_mesh.z.size()), tmp);
	_rho.resize(int(_mesh.z.size()), tmp);
	_p0.resize(int(_mesh.z.size()), tmp);
	_h0.resize(int(_mesh.z.size()), tmp);
	_temp.resize(int(_mesh.z.size()), tmp);
	_rCt.resize(int(_mesh.z.size()), tmp);

	std::fill(tmp.begin(), tmp.end(), 0.);
	_RHS.resize(int(_mesh.z.size()), tmp);
	_entropy.resize(int(_mesh.z.size()), tmp);

	for (int i = 0; i < _nStation; i++)
	{
		_mesh.z[i] = 5.*_chord / (_nStation - 1)*i;
		for (int j = 0; j < _nStream; j++){
			_mesh.r[j] = (_rShd - _rHub) / (_nStream - 1)*j + _rHub;
			_psi[i][j] = (pow(_mesh.r[j], 2.) - pow(_rHub, 2.)) /
				(pow(_rShd, 2.) - pow(_rHub, 2.));
			_rho[i][j] = _rho0In;
			_h0[i][j] = _cp*_T0In;
			_p0[i][j] = _rho0In*_gasConst*_T0In;
			if (i <= _indexLE)
				_rCt[i][j] = _rCtIn;
			if (i >= _indexTE)
				_rCt[i][j] = _rCtOut;
			else
				_rCt[i][j] = (_rCtOut - _rCtIn)*(i - _indexLE) /
				(_indexTE - _indexLE) + _rCtIn;
		}
	}
#if _DEBUG
	std::ofstream myFile;
	myFile.open("./REE/output/checkInitialPsi.csv");
	if (myFile.is_open()) {
		for (int j = _nStream-1; j >= 0; j--) {
			for (int i = 0; i < _nStation; i++) {
				myFile << _psi[i][j] << ",";
			}
			myFile << "\n";
		}
	}
	else {
		std::cout << "Could not open file in initializeField()" << std::endl;
		return false;
	}
#endif 
	return true;
}

bool REESolver::computeVelocityField()
{
	for (int i = 0; i < _nStation; i++) {
		for (int j = 0; j < _nStream; j++){
			if (j == _nStream - 1) {
				_Cz[i][j] = _mfr / (2.*M_PI*_rho[i][j] * _mesh.r[j])*
					(_psi[i][j] - _psi[i][j - 1]) /
					(_mesh.r[j] - _mesh.r[j - 1]);
			}
			else if (j == 0) {
				_Cz[i][j] = _mfr / (2.*M_PI*_rho[i][j] * _mesh.r[j])*
					(_psi[i][j+1] - _psi[i][j]) /
					(_mesh.r[j+1] - _mesh.r[j]);
			}
			else {
				_Cz[i][j] = _mfr / (2.*M_PI*_rho[i][j] * _mesh.r[j])*
					(_psi[i][j + 1] - _psi[i][j - 1]) /
					(_mesh.r[j + 1] - _mesh.r[j - 1]);
			}
			if (i == _nStation - 1) {
				_Cr[i][j] = -_mfr / (2.*M_PI*_rho[i][j] * _mesh.r[j])*
					(_psi[i][j] - _psi[i - 1][j]) /
					(_mesh.z[i] - _mesh.z[i - 1]);
			}
			else if (i == 0) {
				_Cr[i][j] = -_mfr / (2.*M_PI*_rho[i][j] * _mesh.r[j])*
					(_psi[i + 1][j] - _psi[i][j]) /
					(_mesh.z[i + 1] - _mesh.z[i]);
			}
			else {
				_Cr[i][j] = -_mfr / (2.*M_PI*_rho[i][j] * _mesh.r[j])*
					(_psi[i+1][j] - _psi[i-1][j]) /
					(_mesh.z[i + 1] - _mesh.z[i - 1]);
			}
		}
	}
	return true;
}

bool REESolver::computeThermoVariables()
{
	_convCritRho = 0.; 
	/*!****************************************************************
	!TRACING OF THE THERMODYNAMIC VARIABLES STATION - BY - STATION
	!*****************************************************************/
	double cSquare(-1.); double  hStatic(-1.);
	double pStatic(-1.);
	double psiDwn(-1.); double psiUp(-1.);
	double psiTarget(-1.); double slope(-1.);
	double omega(0.);

	double RCU1(-1.); double h01(-1.);
	double p01(-1.); double RAD1(-1.);

	double DENS1(-1.); double CZ1(-1.);
	double CR1(-1.); double entropy1(-1.);
	double C1SQ(-1.); double HSTAT1(-1.); 
	double HOR1(-1.); double POR1(-1.);
	double PSTAT1(-1.);
	double ROTALP1(-1.);

	double h0rel2(-1.); double POR2IDL(-1.);
	double p0rel2(-1.);

	double gammaRatio(_gamma / (_gamma - 1.));

	//double ERRDENS(0.);
	double deltaRho(0.);

	int indexLeft(0);
	int indexStart(0);

	double omegaLoss(-1.); double pLoss(-1.);
	/*-------------------------------
	!UPDATE THE DENSITY ON THE INLET
	!--------------------------------*/
	for (int j = 0; j < _nStream; j++) {
		cSquare = pow(_rCt[0][j] / _mesh.r[j], 2.) + pow(_Cz[0][j], 2.) +
			pow(_Cr[0][j], 2.);
		hStatic = _cp*_T0In - 0.5*cSquare;
		pStatic = _p0In * pow(hStatic / _cp / _T0In, gammaRatio);
		_rho[0][j] = pStatic / (_gasConst*hStatic / _cp);
	}
	/*!
	!SWEEP THE PLANES FROM PLANE 2 TO NSTATN
	!-------------------------------------- -
	!*/
	for (int i = 1; i<_nStation; i++) {
		/*!TRACING EACH STREAMLINE BACK TO THE PREVIOUS STATION
		!FOR A DUCT(REFERENCE PLANE)
		!OR TO THE LEADING EDGE IN CASE OF A BLADE ROW
		!--------------------------------------------------------*/
		indexLeft = i - 1;
		if ((i>_indexLE) && (i <= _indexTE)) {
			indexLeft = _indexLE;
		}

		indexStart = 0;
		psiDwn = _psi[indexLeft][indexStart];
		psiUp = _psi[indexLeft][indexStart + 1];
		/*!
		!SWEEP FROM HUB TO SHROUD AT EACH REFERENCE PLANE
		!TO LOCATE ORIGIN OF STREAMLINE
		!------------------------------------------------*/
		for (int j = 0; j<_nStream; j++) {
			psiTarget = _psi[i][j];

			while (psiTarget>psiUp || psiTarget < psiDwn) {
				indexStart++;
				psiDwn = _psi[indexLeft][indexStart];
				psiUp = _psi[indexLeft][indexStart + 1];
				if (indexStart == _nStream)
					std::cout << "Could not trace back streamline from location " << i << std::endl;
			}
			/*!Origin of streamline succesfully located
			!----------------------------------*/
			slope = (psiTarget - psiDwn) / (psiUp - psiDwn);
			//!
			/*!DETERMINE WHETHER THIS A ROTATING BLADE
			!THEN CALCULATE ROTHALPY AT PREVIOUS STATION
			!IF NOT CALCULATE TOTAL ENTHALPY
			!------------------------------------------ -*/
			if ((i>_indexLE) && (i <= _indexTE)) {
				omega = _omega;
			}
			/*!
			!QUANTITIES AT "1" NEEDED FOR CONSERVATION PRINCIPLE
			!-------------------------------------------------- -
			!*/
			RCU1 = slope * (_rCt[indexLeft][indexStart + 1] - _rCt[indexLeft][indexStart])      
				+_rCt[indexLeft][indexStart];
			h01 = slope * (_h0[indexLeft][indexStart + 1] - _h0[indexLeft][indexStart])
				+ _h0[indexLeft][indexStart];
			p01 = slope * (_p0[indexLeft][indexStart + 1] - _p0[indexLeft][indexStart])
				+ _p0[indexLeft][indexStart];
			RAD1 = slope * (_mesh.r[indexStart + 1] - _mesh.r[indexStart])
				+ _mesh.r[indexStart];

			/*!QUANTITIES AT "1" NEEDED FOR LOSS CALCULATION
			!-------------------------------------------- -*/
			DENS1 = slope*(_rho[indexLeft][indexStart + 1] - _rho[indexLeft][indexStart])
				+ _rho[indexLeft][indexStart];
			CZ1 = slope*(_Cz[indexLeft][ indexStart + 1] - _Cz[indexLeft][ indexStart])
				+ _Cz[indexLeft][ indexStart];
			CR1 = slope*(_Cr[indexLeft][indexStart + 1] - _Cr[indexLeft][indexStart])
				+ _Cr[indexLeft][indexStart];
			entropy1 = slope*(_entropy[indexLeft][indexStart + 1] - _entropy[indexLeft][indexStart])
				+ _entropy[indexLeft][indexStart];
			C1SQ = pow(CZ1, 2.) + pow(CR1, 2.) + pow(RCU1 / RAD1, 2.);
			HSTAT1 = h01 - C1SQ / 2.;
			PSTAT1 = p01 * pow(HSTAT1 / h01, gammaRatio);
			/*!
			!ROTATING AND NON - ROTATING QUANTITIES AT REFERENCE STATION
			!-------------------------------------------------------- -*/
			ROTALP1 = h01 - omega * RCU1;
			h0rel2 = ROTALP1 + pow(omega * _mesh.r[j],2.) / 2.;
			HOR1 = ROTALP1 + pow(omega * RAD1,2.) / 2.;
			POR1 = p01 * pow(HOR1 / h01, gammaRatio);
			POR2IDL = POR1 *pow(h0rel2 / HOR1, gammaRatio);

			if ((i > _indexLE) && (i <= _indexTE)) {
				omegaLoss = _lossesTE*double(i - _indexLE) / double(_indexTE - _indexLE);
				pLoss = omegaLoss * (POR1 - PSTAT1);
			} else {
				omegaLoss = 0.;
				pLoss = 0.;
			}

			p0rel2 = POR2IDL - pLoss;
			
			/*!
			!PROJECTS 1 AND 2 : ROTOR WITH "RCU" SPECIFIED
			!------------------------------------------------------ -
			!*/
			_rCt[i][j] = RCU1;
			if (i == _indexLE) _rCt[i][j] = _rCtIn;
			if (i == _indexTE) _rCt[i][j] = _rCtOut;
			//IF(I.EQ.NLE + 1) RCU(I, J) = 117.85
			//IF(I.EQ.NLE + 2) RCU(I, J) = 157.1
			//IF(I.EQ.NLE + 3) RCU(I, J) = 196.35
			//IF(I.EQ.NTE) RCU(I, J) = 235.6
			_h0[i][j] = h01 + omega * (_rCt[i][j] - RCU1);
			_p0[i][j] = p0rel2 * pow(_h0[i][j] / h0rel2,gammaRatio);

			/*!PROJECT 3 : ROTOR WITH "PO2/PO1" SPECIFIED
			!------------------------------------------------------ -
			!
			!PTOTAL(I, J) = p01 * SPECIFIED_PRESSURE_RATIO
			!HTOTAL(I, J) = h0rel2 * (PTOTAL(I, J) / p0rel2)**gammaRatio
			!RCU(I, J) = RCU1 + (HTOTAL(I, J) - h01) / omega

			/*------------------------
			!COMMON CALCULATION BLOCK
			!-------------------------*/
			double Ct = _rCt[i][j] / _mesh.r[j];
			double Wt = Ct - omega * _mesh.r[j];

			double V2SQ = pow(Wt,2.) + pow(_Cz[i][j],2.) + pow(_Cr[i][j],2.);
			double C2SQ = pow(Ct, 2.) + pow(_Cz[i][j], 2.) + pow(_Cr[i][j], 2.);
			hStatic = h0rel2 - V2SQ / 2.0;
			_h0[i][j] = hStatic + C2SQ / 2.0;
			_p0[i][j] = p0rel2 * pow(_h0[i][j] / h0rel2,gammaRatio);
			pStatic = _p0[i][j] * pow(hStatic / _h0[i][j],gammaRatio);
			_temp[i][j] = hStatic / _cp;
			double oldRho = _rho[i][j];
			_rho[i][j] = pStatic / (_gasConst*_temp[i][j]);
			deltaRho = abs(oldRho - _rho[i][j]);
			_convCritRho = std::max(deltaRho, _convCritRho);
			_entropy[i][j] = _cp * log(_h0[i][j] / h01) -
				_gasConst * log(_p0[i][j] / p01)
				+ entropy1;

		}
	}

#if _DEBUG
	std::ofstream myFile;
	myFile.open("./REE/output/checkRCtDist.csv");
	if (myFile.is_open()) {
		for (int j = _nStream - 1; j >= 0; j--) {
			for (int i = 0; i < _nStation; i++) {
				myFile << _rCt[i][j] << ",";
			}
			myFile << "\n";
		}
		myFile.close();
	}
	else {
		std::cout << "Could not open file in computeThermoVariables()" << std::endl;
		return false;
	}
#endif 
	return true;
}

bool REESolver::computeRHS()
{
	double dr(-1.);
	double oldRHS(-1.);
	double deltaRHS(-1.);
	_convCritRHS = 0.;
	for (int i = 0; i < _nStation; i++) {
		for (int j = 1; j < _nStream - 1; j++) {
			oldRHS = _RHS[i][j];
			dr = _mesh.r[j + 1] - _mesh.r[j - 1];
			_RHS[i][j] = -2.*M_PI / (dr*_mfr*_Cz[i][j])*
				(_rCt[i][j] / pow(_mesh.r[j], 2.)*(_rCt[i][j + 1] - _rCt[i][j - 1]) +
				_temp[i][j] * (_entropy[i][j + 1] - _entropy[i][j - 1]) -
				(_h0[i][j + 1] - _h0[i][j - 1]));
			deltaRHS = abs(_RHS[i][j] - oldRHS);
			_convCritRHS = std::max(deltaRHS, _convCritRHS);
		}
	}
	return true;
}

bool REESolver::updateStreamField()
{
	// Update psi function (stream function) using a "point relaxation technique"
	std::vector<std::vector<double> > A;
	std::vector<std::vector<double> > B;

	std::vector<double> tmp(int(_mesh.r.size()), -1.);
	A.resize(int(_mesh.z.size()), tmp);
	B.resize(int(_mesh.z.size()), tmp);

	double rhorjPlus(-1.);
	double rhorjMinus(-1.);
	double rhoriPlus(-1.);
	double rhoriMinus(-1.);
	double dz(-1.);
	double dr(-1.);
	double oldPsi(-1.);
	double deltaPsi(-1.);
	int count(0);
	int iPlus(0);
	int iMinus(0);

	while (_convCritPsi > _tolPsi || count < 10) {
		_convCritPsi = 0.;
		for (int i = 0; i < _nStation; i++) {
			for (int j = 1; j < _nStream-1; j++) {
				if (i == 0) {
					iMinus = i;
				}
				else {
					iMinus = i - 1;
				}
				if (i == _nStation - 1) {
					iPlus = i;
				}
				else {
					iPlus = i + 1;
				}
				rhorjPlus = 0.5*(_mesh.r[j + 1] * _rho[i][j + 1] + _mesh.r[j] * _rho[i][j]);
				rhorjMinus = 0.5*(_mesh.r[j] * _rho[i][j] + _mesh.r[j - 1] * _rho[i][j - 1]);
				rhoriPlus = 0.5*(_mesh.r[j] * _rho[iPlus][j] + _mesh.r[j] * _rho[i][j]);
				rhoriMinus = 0.5*(_mesh.r[j] * _rho[i][j] + _mesh.r[j] * _rho[iMinus][j]);
				dr = 0.5*(_mesh.r[j + 1] + _mesh.r[j]) - 0.5*(_mesh.r[j] + _mesh.r[j - 1]);
				dz = 0.5*(_mesh.z[iPlus] + _mesh.z[i]) - 0.5*(_mesh.z[i] + _mesh.z[iMinus]);
				A[i][j] = pow(1. / rhoriPlus + 1. / rhoriMinus + pow(dz / dr, 2.) / rhorjPlus +
					pow(dz / dr, 2.) / rhorjMinus, -1.);
				if (i == _nStation - 1) // Neumann Boundary Conditions apply i.e. delPsi/delZ = 0
				{
					/*B[i][j] = _psi[iPlus][j] / rhoriPlus + _psi[iPlus][j] / rhoriMinus +
						pow(dz / dr, 2.)*(_psi[i][j + 1] / rhorjPlus + _psi[i][j - 1] / rhorjMinus);*/
					B[i][j] = 2.*_psi[iPlus][j] / _rho[i][j] / _mesh.r[j] +
						pow(dz / dr, 2.)*(_psi[i][j + 1] / rhorjPlus + _psi[i][j - 1] / rhorjMinus);
				}
				else {
					B[i][j] = _psi[iPlus][j] / rhoriPlus + _psi[iMinus][j] / rhoriMinus +
						pow(dz / dr, 2.)*(_psi[i][j + 1] / rhorjPlus + _psi[i][j - 1] / rhorjMinus);
				}
				oldPsi = _psi[i][j];
				_psi[i][j] = A[i][j] * (B[i][j] + dz*dz*_RHS[i][j]);
				deltaPsi = abs(_psi[i][j] - oldPsi);
				_convCritPsi = std::max(deltaPsi, _convCritPsi);
			}
		}
		count++;
		//std::cout << "sweep #" << count << " Psi max. diff = " << _convCritPsi << std::endl;
		if (count > _maxIt) {
			std::cout << "Maximum number of iterations reached in updateStreamField" << std::endl;
			return false;
		}
	}
#if _DEBUG
	std::ofstream myFile;
	myFile.open("./REE/output/checkUPdatedPsi.csv");
	if (myFile.is_open()) {
		for (int j = _nStream - 1; j >= 0; j--) {
			for (int i = 0; i < _nStation; i++) {
				myFile << _psi[i][j] << ",";
			}
			myFile << "\n";
		}
	}
	else {
		std::cout << "Could not open file" << std::endl;
		return false;
	}
#endif 
	//TODO:
	// implement Gauss-Seidel solver for solving all the points from the mesh in one go
	return true;
}

bool REESolver::writeReport()
{
	double Cm(-1.); 
	std::ofstream myFile;
	myFile.open("./REE/output/REE-Report.csv");
	if (myFile.is_open()) {
		myFile << "The following values are measured at TE\n";
		myFile << "Cz,Cr,beta,rho,r\n";
		for (int j = _nStream - 1; j >= 0; j--) {
			Cm = pow(pow(_Cz[_indexTE][j],2)+pow(_Cr[_indexTE][j],2.),0.5);
			_beta[j] = atan(abs(_mesh.r[j]*_omega-
				_rCt[_indexTE][j]/_mesh.r[j])/Cm);
			myFile << _Cz[_indexTE][j] << "," << _Cr[_indexTE][j] << "," <<
				_beta[j]*180./M_PI << "," << _rho[_indexTE][j] << "," << 
				_mesh.r[j] << "\n"; 
		}
		myFile.close();
	}
	else {
		std::cout << "Could not open file writeReport()" << std::endl;
		return false;
	}
	return true; 
}



