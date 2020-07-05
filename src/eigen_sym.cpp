#include "eigen_sym.h"

#include <cmath>
#include <iostream>

using namespace std;

bool eig_sym(REAL *pMatrix, INTG nDim, REAL *eigVal, REAL *eigVec, REAL dbEps, INTG nJt)
{
    bool    ret = false;

    if (nJt == -1)
        nJt = nDim * nDim * 100;

    for(INTG i=0; i<nDim; i++)
        for (INTG j=0; j<nDim; j++)
            eigVec[i*nDim+j] = (i == j ? 1.0 : 0.0);


    INTG nCount = 0;
    while(true) {
        REAL dbMax = pMatrix[1];
        INTG nRow = 0;
        INTG nCol = 1;
        for (INTG i=0; i<nDim; i++) {
            for (INTG j=0; j<nDim; j++) {
                REAL d = abs(pMatrix[i * nDim + j]);
                if(i != j && d > dbMax) {
                    dbMax = d;
                    nRow  = i;
                    nCol  = j;
                }
            }
        }

        if(dbMax < dbEps) {
            ret = true;
            break;
        }

        if(nCount > nJt) {
            ret = false;
            break;
        }

        nCount++;

        REAL dbApp = pMatrix[nRow * nDim + nRow];
        REAL dbApq = pMatrix[nRow * nDim + nCol];
        REAL dbAqq = pMatrix[nCol * nDim + nCol];

        // calculate angle
        REAL dbAngle        = 0.5*atan2(-2*dbApq,dbAqq-dbApp);
        REAL dbSinTheta     = sin(dbAngle);
        REAL dbCosTheta     = cos(dbAngle);
        REAL dbSin2Theta    = sin(2*dbAngle);
        REAL dbCos2Theta    = cos(2*dbAngle);

        pMatrix[nRow*nDim+nRow] = dbApp*dbCosTheta*dbCosTheta +
                                  dbAqq*dbSinTheta*dbSinTheta +
                                  2*dbApq*dbCosTheta*dbSinTheta;
        pMatrix[nCol*nDim+nCol] = dbApp*dbSinTheta*dbSinTheta +
                                  dbAqq*dbCosTheta*dbCosTheta -
                                  2*dbApq*dbCosTheta*dbSinTheta;
        pMatrix[nRow*nDim+nCol] = 0.5*(dbAqq-dbApp)*dbSin2Theta + dbApq*dbCos2Theta;
        pMatrix[nCol*nDim+nRow] = pMatrix[nRow*nDim+nCol];

        for(INTG i = 0; i < nDim; i ++) {
            if((i!=nCol) && (i!=nRow)) {
                INTG u = i*nDim + nRow;
                INTG w = i*nDim + nCol;
                dbMax = pMatrix[u];
                pMatrix[u]= pMatrix[w]*dbSinTheta + dbMax*dbCosTheta;
                pMatrix[w]= pMatrix[w]*dbCosTheta - dbMax*dbSinTheta;
            }
        }

        for (INTG j = 0; j < nDim; j ++) {
            if((j!=nCol) && (j!=nRow)) {
                INTG u = nRow*nDim + j;
                INTG w = nCol*nDim + j;
                dbMax = pMatrix[u];
                pMatrix[u]= pMatrix[w]*dbSinTheta + dbMax*dbCosTheta;
                pMatrix[w]= pMatrix[w]*dbCosTheta - dbMax*dbSinTheta;
            }
        }

        for(int i = 0; i < nDim; i ++)
        {
            INTG u = i*nDim + nRow;
            INTG w = i*nDim + nCol;
            dbMax = eigVec[u];
            eigVec[u] = eigVec[w]*dbSinTheta + dbMax*dbCosTheta;
            eigVec[w] = eigVec[w]*dbCosTheta - dbMax*dbSinTheta;
        }

    } // end of while

    for (INTG i=0; i<nDim; i++)
        eigVal[i] = pMatrix[i * nDim + i];

    return ret;
}
