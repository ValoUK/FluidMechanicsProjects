
/*!****************************************************************
!TRACING OF THE THERMODYNAMIC VARIABLES STATION - BY - STATION
!*****************************************************************/
double CSQ(-1.); double  HSTATIC(-1.);
double PSTATIC(-1.);
double psiDwn(-1.); double psiUp(-1.);
double psiTarget(-1.); double slope(-1.);
double ROTATE(-1.);

double RCU1(-1.); double HTOTAL1(-1.);
double PTOTAL1(-1.); double RAD1(-1.);

double DENS1(-1.); double CZ1(-1.);
double CR1(-1.); double ENTROP1(-1.);
double C1SQ(-1.); double HSTAT1(-1.);
double PSTAT1(-1.);

double GAMAM(_gamma / (_gamma - 1.));
/*!COMMON / NUMBER / H,
!
!UPDATE THE DENSITY ON THE INLET
!------------------------------ -*/
for (int j = 0; j < _nStation; j++) {
	CSQ = pow(_rCt[0][j] / _mesh.r[j], 2.) + _Cz[0][j] * _Cz[0][j] +
		_Cr[0][j] * _Cr[0][j]
		HSTATIC = _h0[0][j] - 0.5*CSQ;
	PSTATIC = _p0[0][j] * pow(HSTATIC / _h0[0][j], _gamma*_gasConst / (_gamma - 1.));
	_rho[0][j] = PSTATIC / (_gasConst*HSTATIC / _cp)
}
/*!
!SWEEP THE PLANES FROM PLANE 2 TO NSTATN
!-------------------------------------- -
!*/
double ERRDENS = 0.0
double CHANGE = 0.0

//!-- -

for (int i = 1; i<_nStation; i++) {//DO I = 2, NSTATN
	/*!TRACING EACH STREAMLINE BACK TO THE PREVIOUS STATION
	!FOR A DUCT(REFERENCE PLANE)
	!OR TO THE LEADING EDGE IN CASE OF A BLADE ROW
	!--------------------------------------------------------*/
	_indexLeft = I - 1
		if (i>_indexLE && i <= _indexTE) {//IF(I.GT.NLE.AND.I.LE.NTE) THEN
			_indexLeft = _indexLE
		}

	_indexStart = 0;
	psiDwn = PSI(_indedLeft, _indexSTart)
		psiUp = PSI(_indexLeft, _indexSTart + 1)
		/*!
		!SWEEP FROM HUB TO SHROUD AT EACH REFERENCE PLANE
		!TO LOCATE ORIGIN OF STREAMLINE
		!------------------------------------------------*/
		for (int j = 0; j<_nStream; j++) {//DO  J = 1, NSTRM
			//!
			psiTarget = _psi[i][j];
			//!
			while (psiTarget>psiUp || psiTarget < psiDwn) {
				_indexStart++;
				psiDwn = _psi[_indexLeft][_indexStart];
				psiUp = _psi[_indexLeft][_indexStart + 1];
				if (_indexStart == _nStream)
					std::cout << "Could not trace streamline" << std::endl;
			}
			/*!STREAMLINE ORIGIN HAS BEEN LOCATED
			!----------------------------------*/
			slope = (psiTarget - psiDwn) / (psiUp - psiDwn);
			//!
			/*!DETERMINE WHETHER THIS A ROTATING BLADE
			!THEN CALCULATE ROTHALPY AT PREVIOUS STATION
			!IF NOT CALCULATE TOTAL ENTHALPY
			!------------------------------------------ -*/
			ROTATE = 0.0;
			if (i>_indexLE && i <= _indexTE) {
				ROTATE = _omega;
			}
			/*!
			!QUANTITIES AT "1" NEEDED FOR CONSERVATION PRINCIPLE
			!-------------------------------------------------- -
			!*/
			RCU1 = slope * (RCU(_indexLeft, NSTART + 1) - RCU(_indexLeft, NSTART))      &
				+RCU(_indexLeft, NSTART)
				HTOTAL1 = slope * (HTOTAL(_indexLeft, NSTART + 1) - HTOTAL(_indexLeft, NSTART))    &
				+HTOTAL(_indexLeft, NSTART)
				PTOTAL1 = slope * (PTOTAL(_indexLeft, NSTART + 1) - PTOTAL(_indexLeft, NSTART))    &
				+PTOTAL(_indexLeft, NSTART)
				RAD1 = slope * (RADIUS(_indexLeft, NSTART + 1) - RADIUS(_indexLeft, NSTART))       &
				+RADIUS(_indexLeft, NSTART)
				//!

				//!---- -

				/*!QUANTITIES AT "1" NEEDED FOR LOSS CALCULATION
				!-------------------------------------------- -*/
				DENS1 = slope*(DENSITY(_indexLeft, NSTART + 1) - DENSITY(_indexLeft, NSTART))
				+ DENSITY(_indexLeft, NSTART);
			CZ1 = slope*(CZ(_indexLeft, NSTART + 1) - CZ(_indexLeft, NSTART))
				+ CZ(_indexLeft, NSTART);
			CR1 = slope*(CR(_indexLeft, NSTART + 1) - CR(_indexLeft, NSTART))
				+ CR(_indexLeft, NSTART);
			ENTROP1 = slope*(ENTROPY(_indexLeft, NSTART + 1) - ENTROPY(_indexLeft, NSTART))
				+ ENTROPY(_indexLeft, NSTART);
			C1SQ = pow(CZ1, 2.) + pow(CR1, 2.) + pow(RCU1 / RAD1, 2.);
			HSTAT1 = HTOTAL1 - C1SQ / 2.;
			PSTAT1 = PTOTAL1 * pow(HSTAT1 / HTOTAL1, GAMAM);
			/*!
			!ROTATING AND NON - ROTATING QUANTITIES AT REFERENCE STATION
			!-------------------------------------------------------- -*/
			ROTALP1 = HTOTAL1 - ROTATE * RCU1;
			HOR2 = ROTALP1 + (ROTATE * RADIUS(I, J))**2 / TWO;
			HOR1 = ROTALP1 + (ROTATE * RAD1) **2 / TWO;
			POR1 = PTOTAL1 * pow(HOR1 / HTOTAL1, GAMAM);
			POR2IDL = POR1 *pow(HOR2 / HOR1, GAMAM);

			IF(I.GT.NLE.AND.I.LE.NTE) THEN
				omegaLoss = 0.03*FLOAT(I - NLE) / FLOAT(NTE - NLE)
				pLoss = omegaLoss * (POR1 - PSTAT1)
				ELSE
				omegaLoss = 0.0
				pLoss = 0.0
				END IF

				POR2 = POR2IDL - pLoss
				/*!
				!PROJECTS 1 AND 2 : ROTOR WITH "RCU" SPECIFIED
				!------------------------------------------------------ -
				!*/
				/*RCU(I, J) = RCU1
				IF(I.EQ.NLE + 1) RCU(I, J) = 117.85
				IF(I.EQ.NLE + 2) RCU(I, J) = 157.1
				IF(I.EQ.NLE + 3) RCU(I, J) = 196.35
				IF(I.EQ.NTE) RCU(I, J) = 235.6
				HTOTAL(I, J) = HTOTAL1 + ROTATE * (RCU(I, J) - RCU1)
				PTOTAL(I, J) = POR2 * (HTOTAL(I, J) / HOR2)**GAMAM*/

				/*!PROJECT 3 : ROTOR WITH "PO2/PO1" SPECIFIED
				!------------------------------------------------------ -
				!
				!PTOTAL(I, J) = PTOTAL1 * SPECIFIED_PRESSURE_RATIO
				!HTOTAL(I, J) = HOR2 * (PTOTAL(I, J) / POR2)**GAMAM
				!RCU(I, J) = RCU1 + (HTOTAL(I, J) - HTOTAL1) / ROTATE

				/*------------------------
				!COMMON CALCULATION BLOCK
				!-------------------------*/
				/*CU = RCU(I, J) / RADIUS(I, J)
				VU = CU - ROTATE * RADIUS(I, J)*/

				/*V2SQ = VU**2 + CZ(I, J)**2 + CR(I, J)**2
				C2SQ = CU**2 + CZ(I, J)**2 + CR(I, J)**2
				HSTATIC = HOR2 - V2SQ / 2.0
				HTOTAL(I, J) = HSTATIC + C2SQ / 2.0
				PTOTAL(I, J) = POR2 * (HTOTAL(I, J) / HOR2)**GAMAM
				PSTATIC = PTOTAL(I, J) * (HSTATIC / HTOTAL(I, J))**GAMAM
				DENSOLD = DENSITY(I, J)
				DENSITY(I, J) = PSTATIC / (RGAS*HSTATIC / CP)
				CHANGE = ABS(DENSOLD - DENSITY(I, J))
				ERRDENS = AMAX1(CHANGE, ERRDENS)
				ENTROPY(I, J) = CP * LOG(HTOTAL(I, J) / HTOTAL1) - &
				RGAS * LOG(PTOTAL(I, J) / PTOTAL1)                 &
				+ENTROP1*/

		}
}
	/*!
	!ASSEMBLE THE RHS FOR UNKNOWN NODES
	!I.E.ALL NODES EXCEPT HUB / SHROUD
	!----------------------------------
	!*/

	/*ERRRHS = 0.0
	CHANGE = 0.0

	DO  I = 1, NSTATN
		DO  J = 2, NSTRM - 1

			RHSOLD = RHS(I, J)
			CU = RCU(I, J) / RADIUS(I, J)
			TSTATIC = (HTOTAL(I, J) - (CZ(I, J)**2 + CR(I, J)**2 + CU**2) / 2.0) / CP
			RHS(I, J) = -1.0 / CZ(I, J) * 2.0*PI / XMASS *            &
				(CU / RADIUS(I, J) *(RCU(I, J + 1) - RCU(I, J - 1)) / TWOH  &
				-(HTOTAL(I, J + 1) - HTOTAL(I, J - 1)) / TWOH             &
				+TSTATIC*(ENTROPY(I, J + 1) - ENTROPY(I, J - 1)) / TWOH) &
			CHANGE = ABS(RHS(I, J) - RHSOLD)
			ERRRHS = AMAX1(ERRRHS, CHANGE)

		END DO
	END DO*/

	//999 FORMAT(10X, 'PSI VALUE = ', F10.3, 'STATION I =', I5, 2X, &
	//'STREAMLINE J =', I5, / , 'PSI UP = ', F10.3, 'PSI DOWN = ', F10.3)

	//!RETURN