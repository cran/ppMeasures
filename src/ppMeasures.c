#include <R.h>
#include <Rmath.h>


void fixMatching(double *AA, int *n1, int *n2, int *xtoy,
				 int *ytox, int *foundCycle, double *eps)
{ /* we look for cycles and eliminate them */
	int i=2, tI=-1, j=0, tJ=-1, holdToI;
	double theBest=1, adds=5, subs=1;
	
	for(i=0; i<*n1; i++){
		if(xtoy[i] > -1){ 
			for(j=i+1; j<*n1; j++){
				if(xtoy[j] > -1){
					/* try the switch */
					subs = AA[i* *n2 + xtoy[i]] + AA[j* *n2 + xtoy[j]];
					adds = AA[i* *n2 + xtoy[j]] + AA[j* *n2 + xtoy[i]];
					if(adds - subs < theBest - *eps){
						theBest = adds - subs;
						tI = i;
						tJ = j;
					}
				}
			}
		}
	}
	if(theBest < 0){
		holdToI = xtoy[tI];
		xtoy[tI] = xtoy[tJ];
		xtoy[tJ] = holdToI;
		
		holdToI = ytox[xtoy[tI]];
		ytox[xtoy[tI]] = ytox[xtoy[tJ]];
		ytox[xtoy[tJ]] = holdToI;
		
		*foundCycle=1;
	}
}

void stdist(double *p1,  double *p2,  int *n1,    int *n2,   int *dim,
			double *pm,  double *pa,  double *pd, double *eps, int *xtoy,
			int *ytox, int *maxB, double *cost, int *nB, int *euclid,
			double *lossOrder, int *kernProto)
{
	int    d, p, i, I, j, jj, J, JJ, k, K, tI, tJ, tK,
		stillMatching=1, updates, foundCycle;
	double P1[*n1][*dim], P2[*n2][*dim],
		A[*n1][*n2], An[*n1][*n2], Ao[*n1][*n2], AA[*n1* *n2],
		temp, temp1, dMin;
	
	for(p=0; p<*n1; p++){
		for(d=0; d<*dim; d++){
			P1[p][d] = p1[*n1 * d + p];
		}
		xtoy[p] = -1;
	}
	for(p=0; p<*n2; p++){
		for(d=0; d<*dim; d++){
			P2[p][d] = p2[*n2 * d + p];
		}
		ytox[p] = -1;
	}
	K=0;
	for(i=0; i<*n1; i++){
		for(j=0; j<*n2; j++){
			temp = 0;
			temp1 = 0;
			for(k=0; k<*dim; k++){
				if(euclid[k] == 0){
					temp += pm[k] * pow(fabs(P1[i][k] - P2[j][k]), *lossOrder);
				} else {
					//temp1 += pow(pm[k], 2/ *lossOrder)*(P1[i][k]-P2[j][k])*(P1[i][k]-P2[j][k]);
					temp1 += pm[k]*pm[k]*(P1[i][k]-P2[j][k])*(P1[i][k]-P2[j][k]);
				}
			}
			A[i][j] = temp + pow(sqrt(temp1), *lossOrder);
			AA[i* *n2 + j] = A[i][j];
		}
	}
	if(*n1 < 2 || *n2 < 2){
		*cost = *n1 * *pd + *n2 * *pa;
		if(*n1 == 1 && *n2 > 0){
			J = 0;
			for(j=1; j<*n2; j++){
				if(A[0][j] < A[0][J]){
					J = j;
				}
			}
			if(A[0][J] < *pa+*pd){
				xtoy[0] = J;
				ytox[J] = 0;
				*cost = (*n2-1) * *pa + A[0][J];
			}
		} else if(*n2 == 1 && *n1 > 0){
			I = 0;
			for(i=1; i<*n1; i++){
				if(A[i][0] < A[I][0]){
					I = i;
				}
			}
			if(A[I][0] < *pa+*pd){
				xtoy[I] = 0;
				ytox[0] = I;
				*cost = (*n1-1) * *pd + A[I][0];
			}
		}
	} else {
		while(stillMatching == 1 && K < *n1 && K < *n2){
			if(K > 0){
				tK = 0;
				tI = 0;
				tJ = 0;
				dMin = *pa + *pd + 1; 
				for(i=0; i<*n1; i++){
					if(xtoy[i] == -1){
						for(j=0; j<*n2; j++){
							An[i][j] = A[i][j];
							if(A[i][j] < dMin && ytox[j] == -1){
								dMin = A[i][j];
								tI = i;
								tJ = j;
							}
						}
					}
				}
				updates = 1;
				k       = 1;
				/* An[i][j] (at iteration k) represents the cost of matching
				 i and j while breaking up to k current matchings */
				while(k < K+1 && k < *maxB+1 && updates > 0){
					for(i=0;i<*n1;i++){
						for(j=0;j<*n2;j++){
							Ao[i][j] = An[i][j];
						}
					}
					updates = 0;
					for(i=0;i<*n1;i++){
						for(j=0;j<*n2;j++){
							if(xtoy[i] == -1){
								An[i][j] = Ao[i][j];
								JJ=-1;
								for(jj=0; jj<*n2; jj++){
									if(ytox[jj] > -1){
										temp = Ao[i][jj] - A[ytox[jj]][jj] + A[ytox[jj]][j];
										if(temp < An[i][j] - *eps){
											An[i][j] = temp;
											JJ = jj;
										}
									}
								}
								if(JJ > -1){
									if(ytox[j] == -1){
										if(An[i][j] < dMin - *eps){
											dMin = An[i][j];
											tK = k;
											tI = i;
											tJ = j;
										}
									} else {
										if(An[i][j] < dMin){
											updates = 1;
										}
									}
								}
							} else {
								An[i][j] = *pa + *pd + 1;
							}
						}
					}
					k++;
				}
				if(dMin < *pa + *pd){
					xtoy[tI] = tJ;
					ytox[tJ] = tI;
					if(tK > 0 || *maxB == 0){
						foundCycle = 1;
						if(*maxB == 0 && *kernProto == 1){
							foundCycle = 0;
						}
						while(foundCycle == 1){
							foundCycle = 0; 
							fixMatching(AA, n1, n2, xtoy, ytox, &foundCycle, eps);
						}
					}
				} else {
					stillMatching = 0;
					K--;
				}
			} else {
				dMin = A[0][0];
				I = 0;
				J = 0;
				for(i=0;i<*n1;i++){
					for(j=0;j<*n2;j++){
						if(A[i][j] < dMin){
							I = i;
							J = j;
							dMin = A[i][j];
						}
					}
				}
				if(dMin < *pa + *pd){
					xtoy[I] = J;
					ytox[J] = I;
				} else {
					stillMatching = 0;
					K--;
				}
			}
			K++;
		}
		*cost = *n1 * *pd + *n2 * *pa;
		for(i=0; i<*n1; i++){
			if(xtoy[i] > -1){
				*cost += A[i][xtoy[i]] - *pa - *pd;
			}
		}
	}
}



void getPTC(double *pts, int *npts, int *dim,
	    double *PTC, int *nPTC, int *PTCstart, double *posEps)
{
	int d, p, k, accept;
	*nPTC=0;
	for(d=0; d<*dim; d++){
		PTCstart[d] = *nPTC;
		PTC[*nPTC] = pts[*npts * d];
		*nPTC += 1;
		for(p=1; p<*npts; p++){
			PTC[*nPTC] = pts[*npts * d + p];
			accept=1;
			for(k=PTCstart[d]; k<*nPTC; k++){
				if(fabs(PTC[k] - PTC[*nPTC]) < posEps[d]){
					accept=0;
				}
			}
			if(accept == 1){
			  *nPTC += 1;
			}
		}
	}
}

void margKernWMatch(double *pts, int *npts, int *key,
			  int *match, int *keyMax, int *dim,
			  double *PTC, int *nPTC, int *PTCstart,
			  double *pm,  double *pa,  double *pd,
			  double *kern,
			  double *wts, double *lossOrder)
{
	int    d, k, at, stop;
	double tempKern, theBest[*keyMax];
	for(d=0; d<*dim; d++){
		if(d == *dim - 1){
			stop = *nPTC;
		} else {
			stop = PTCstart[d+1];
		}
		for(at=PTCstart[d]; at<stop; at++){ 
			for(k=0; k<*keyMax; k++){ // theBest[k] is the best kern for key k
				theBest[k] = *pd;
			}
			for(k=0; k<*npts; k++){
				tempKern = pm[d]*fabs(PTC[at] - pts[*npts*d + k]) - *pa;
				if(tempKern < theBest[key[k]-1] && match[k] == -1){
					theBest[key[k]-1] = tempKern; // key starts at 1 and goes to keyMax
				}
			}
			kern[at] = wts[0] * theBest[0];
			for(k=1; k<*keyMax; k++){
				kern[at] += wts[k] * theBest[k];
			}
		}
	}
}

void sample(double *pts, int *npts,
			double *Lbd, double *Ubd,
			int *samp, int *nsamp, int *totSampSize)
{
	int d, p, i, j, ind, N, maxI, alreadySampled;
	double pt, dist, ptMin, max;
	if(*nsamp==0){
		// sample the first point in pts that falls in Lbd and Ubd
		for(p=0; p<*npts; p++){
			if(pts[p] < *Ubd && pts[p] > *Lbd){
				samp[*nsamp] = p;
				*nsamp += 1;
				break;
			}
		}
	} // get first sample, if needed
	for(N=*nsamp; N<*totSampSize; N++){
		max = 0;
		for(p=0; p<*npts; p++){
			alreadySampled = 0;
			for(i=0; i<*nsamp; i++){
				if(p == samp[i]){
					alreadySampled = 1;
				}
			} // check if this pt was already sampled
			if(!alreadySampled){
				pt = pts[p];
				if(pt < *Ubd && pt > *Lbd){
					dist  = fabs(pt - pts[samp[0]]);
					ptMin = dist;
					for(i=1; i<*nsamp; i++){
						dist = fabs(pt - pts[samp[i]]);
						if(dist < ptMin){
							ptMin = dist;
						}
					} // get min distance to nearest selected sample pt
					if(ptMin > max){
						max  = ptMin;
						maxI = p;
					} // found a better point
				}
			}
		}
		// if found a new point, then keep it and continue
		// otherwise stop looking for new sample locations
		if(max > 0){
			samp[*nsamp] = maxI;
			*nsamp += 1;
		} else {
			break;
		}
	}
}

void findLattice(double *Lbd, double *Ubd, int *dim,
				 double *PTC,  int *nPTC, int *PTCstart,
				 double *kern,
				 int *ppd, int *maxppd, double *space,
				 double *L, int *nL)
{
	int d, p, P, Lp, i, ptsToGet, NL, same, start, stop,
	nLi[*dim], Li[*maxppd][*dim], tooClose,
	NPTS, SAMP[*maxppd], NSAMP, TOTSAMPSIZE,
	inc[*dim], inc0, go;
	double LBD, UBD, sp[*dim], kMin, kMinIndex, PTS[*nPTC];
	for(d=0; d<*dim; d++){
		sp[d]   = (Ubd[d] - Lbd[d]) * space[d];
		for(p=0; p<*maxppd; p++){
			Li[p][d] = -1;
		}
		nLi[d]   = 0;
	} // initialize sp, Li, and nLi
	for(d=0; d<*dim; d++){
		ptsToGet=0;
		start=PTCstart[d];
		if(d < *dim-1){
			stop = PTCstart[d+1];
		} else {
			stop = *nPTC;
		} // find the stopping place for PTC
		kMinIndex = -1;
		for(p=start; p<stop; p++){
			if(PTC[p] > Lbd[d] && PTC[p] < Ubd[d]){
				ptsToGet++;
				kMin      = kern[p];
				kMinIndex = p;
			}
		} // determine the number of potential locations
		if(ptsToGet <= ppd[d]){ // get all the points (getPTC results in no doubles)
			for(p=start; p<stop; p++){
				if(PTC[p] >= Lbd[d] && PTC[p] <= Ubd[d]){
					Li[nLi[d]][d] = p;
					nLi[d]++;
				}
			}
		} else {
			// kMin & kMinIndex adjusted above to reflect some value >Lbd[d] <Ubd[d]
			for(p=start; p<stop; p++){
				if(kern[p] < kMin && PTC[p] > Lbd[d] && PTC[p] < Ubd[d]){
					kMin      = kern[p];
					kMinIndex = p;
				}
			}
			// note the first minimum
			Li[nLi[d]][d] = kMinIndex;
			nLi[d] += 1;
			//if(kMinIndex < 0){
			//	Rprintf("Finding no minimum in 04findLattice.\n");
			//}
			// get remainder of the points
			for(Lp=1; Lp<ppd[d]; Lp++){
				kMin = 0;
				kMinIndex=-1;
				for(p=start; p<stop; p++){
					if((kern[p] < kMin || kMinIndex == -1)
					   && PTC[p] > Lbd[d] && PTC[p] < Ubd[d]){
						tooClose = 0;
						for(i=0; i<Lp; i++){
							if(fabs(PTC[p] - PTC[Li[i][d]]) < sp[d]){
								tooClose = 1;
							}
						} // verify not too close to other selected pts
						if(!tooClose){
							kMin = kern[p];
							kMinIndex = p;
						} // point becomes a candidate
					}
				} // find the next min, if it is not within sp[d]
				// if there is a pt in region but not within sp of others, add it to Li
				// otherwise, sample from the remaining pts
				if(kMinIndex > -1){
					// update Li & nLi
					Li[nLi[d]][d]  = kMinIndex;
					nLi[d]++;
				} else {
					// sample remaining points
					for(p=0; p<stop-start; p++){
						PTS[p] = PTC[p+start];
					}
					LBD = Lbd[d];
					UBD = Ubd[d];
					NPTS = stop-start;
					for(p=0; p<nLi[d]; p++){
						SAMP[p]=Li[p][d]-start;
					}
					NSAMP = nLi[d];
					TOTSAMPSIZE = ppd[d];
					sample(PTS, &NPTS, &LBD, &UBD, SAMP, &NSAMP, &TOTSAMPSIZE);
					for(p=nLi[d]; p<NSAMP; p++){
						Li[p][d]=SAMP[p]+start;
					} // update Li
					nLi[d] = NSAMP;
					break;
				} 
			}
		}
	}
	// cross product all of the elements of Li
	for(d=0; d<*dim; d++){
		inc[d]=0;
	}
	*nL=0;
	while(inc[0] < nLi[0]){
		for(d=0; d<*dim; d++){
			L[*dim * *nL + d] = PTC[Li[inc[d]][d]];
		}
		*nL += 1;
		go=1;
		for(d=*dim-1; d>-1; d--){
			if(go){
				inc[d] += 1;
				if(inc[d] < nLi[d]){
					go=0;
				} else if(d > 0){
					inc[d]=0;
				}
			}
		}
	}
}



void checkLattice(double *pts, int *npts, int *key,
				  int *keyMax, int *maxObs, int *dim,
				  double *L, int *nL, double *Lcost,
				  double *pm,  double *pa,  double *pd,
				  double *wts, int *maxBranch,
				  double *PT,  int *nPT, double *eps,
				  int *euclid, double *lossOrder, int *useKernOnly,
				  int *kernProto)
{
	int    i, d, p, pp, k, npropPT, npatt, Npatt,
	xtoy[*nPT+1], ytox[*maxObs], nB, theMin, temp;
	double propPT[(*nPT+1)* *dim], patt[*maxObs * *dim],
	cost, Cost;
	npropPT = *nPT+1;
	for(p=0; p<*nPT; p++){
		for(d=0; d<*dim; d++){
			propPT[npropPT*d + p] = PT[p* *dim + d];
		}
	}
	for(pp=0; pp<*nL; pp++){
		for(d=0; d<*dim; d++){
			propPT[npropPT*d + npropPT-1] = L[*dim*pp + d];
		}
		cost=0;
		for(k=0; k<*keyMax; k++){
			// get the # of patterns to check
			Npatt=0;
			for(p=0; p<*npts; p++){
				if(key[p] == k+1){
					Npatt++;
				}
			}
			npatt = 0;
			for(p=0; p<*npts; p++){
				if(key[p] == k+1){
					for(d=0; d<*dim; d++){
						patt[Npatt * d + npatt] = pts[*npts * d + p];
					}
					npatt++;
				}
			}
			// compute costs for each pattern
			if(npatt > 0){
				Cost = 0;
				theMin = npropPT;
				if(npropPT > npatt){
					theMin = npatt;
				}
				nB = 2*(npropPT-1)*(npatt-1);
				for(i=2; i<theMin; i++){
					temp = 2*(npropPT-i)*(npatt-i)*i;
					if(temp > nB){
						nB = temp;
					}
				}
				stdist(propPT, patt, &npropPT, &npatt, dim,
					   pm, pa, pd, eps, xtoy,
					   ytox, maxBranch, &Cost, &nB, euclid,
					   lossOrder, kernProto);
			} else {
				Cost = *pd * npropPT + *pa * Npatt;
			}
			cost += wts[k]*Cost;
		}
		// put costs into Lcost
		Lcost[pp] = cost;
	}
}

void adjustMatch(double *pts, int *npts, int *key,
		   int *match,
		   int *dim, int *keyMax, int *maxObs,
		   double *PT, int *nPT,
		   double *pm, double *pa, double *pd,
		   int *maxBranch, double *wts, double *eps,
		   int *euclid, double *lossOrder,
		   int *kernProto)
{
	int i, k, p, d, N, npatt, Npatt, pattKey[*maxObs],
		xtoy[*nPT], ytox[*maxObs], ind[*maxObs],
		keepPts[*npts], temp, nB, theMin;
	double cost, patt[*maxObs* *dim], Pts[*dim][*maxObs];
	
	// initialize to keep all pts
	for(p=0; p<*npts; p++){
		keepPts[p] = 1;
		match[p]   = -1;
	}
	for(k=0; k<*keyMax; k++){
		// create patt
		Npatt=0;
		for(p=0; p<*npts; p++){
			if(key[p] == k + 1){
				Npatt++;
			}
		}
		npatt=0;
		for(p=0; p<*npts; p++){
			if(key[p] == k + 1){
				for(d=0; d<*dim; d++){
					patt[Npatt*d+npatt] = pts[*npts*d+p];
				}
				ind[npatt] = p;
				npatt++;
			}
		}
		for(p=0; p<*nPT; p++){
			xtoy[p] = -1;
		}
		for(p=0; p<*maxObs; p++){
			ytox[p] = -1;
		}
		
		cost = 0;
		theMin = *nPT;
		if(*nPT > Npatt){
			theMin = Npatt;
		}
		nB = 2*(*nPT-1)*(Npatt-1);
		for(i=2; i<theMin; i++){
			temp = 2*(*nPT-i)*(Npatt-i)*i;
			if(temp > nB){
				nB = temp;
			}
		}
		// match pattern to the PT
		stdist(PT,  patt, nPT, &Npatt, dim,
			   pm, pa, pd, eps, xtoy,
			   ytox, maxBranch, &cost, &nB,
			   euclid, lossOrder, kernProto);
		
		// adjust match to reflect matching
		for(p=0; p<npatt; p++){
			match[ind[p]] = ytox[p];
		}
	}
}

void addPTPt(double *pts, int *npts, int *key,
			 int *match,
			 int *keyMax, int *dim,  int *maxObs,
			 double *PTC, int *nPTC, int *PTCstart,
			 double *pm,  double *pa,  double *pd,
			 double *wts, int *maxBranch,
			 double *PT,  int *nPT, double *eps,
			 int *ppd, int *LMaxElmts, double *space,
			 double *Lbd, double *Ubd, double *newCost,
			 int *euclid, double *lossOrder, int *useKernOnly,
			 int *kernProto)
{
	int d, p, l, LminI, nL, maxppd, stillSearching, stop;
	double Lmin, Lbest[*dim], LBD[*dim], UBD[*dim],
	L[*LMaxElmts* *dim], kern[*nPTC], Lcost[*LMaxElmts],
	lbd, ubd;
	
	// get kernel
	margKernWMatch(pts, npts, key,
			 match, keyMax, dim,
			 PTC, nPTC, PTCstart,
			 pm, pa, pd, kern, wts, lossOrder);
	
	// find the lattice
	maxppd = ppd[0];
	for(d=1; d<*dim; d++){
		if(maxppd < ppd[d]){
			maxppd = ppd[d];
		}
	}
	findLattice(Lbd, Ubd, dim,
				PTC, nPTC, PTCstart, kern, ppd,
				&maxppd, space, L, &nL);
	stillSearching=1;
	stop=0;
	while(stillSearching == 1){
		// check the lattice
		checkLattice(pts, npts, key, keyMax,
					 maxObs, dim, L, &nL, Lcost,
					 pm, pa, pd, wts, maxBranch,
					 PT, nPT, eps, euclid, lossOrder,
					 useKernOnly, kernProto);
		// find the min
		Lmin=Lcost[0];
		LminI=0;
		for(p=1; p<nL; p++){
			if(Lcost[p] < Lmin){
				Lmin=Lcost[p];
				LminI=p;
			}
		}
		for(d=0; d<*dim; d++){
			Lbest[d] = L[*dim *LminI+d];
		}
		// trim bounds: Lbd and Ubd
		for(d=0; d<*dim; d++){
			LBD[d] = Lbd[d];
			for(p=0; p<nL; p++){
				lbd = (Lbest[d] + L[*dim*p+d])/2;
				if(lbd > LBD[d] && L[*dim*p+d] < Lbest[d]){
					LBD[d] = lbd;
				}
			}
			UBD[d] = Ubd[d];
			for(p=0; p<nL; p++){
				ubd = (Lbest[d] + L[*dim*p+d])/2;
				if(ubd < UBD[d] && L[*dim*p+d] > Lbest[d]){
					UBD[d] = ubd;
				}
			}

			if(UBD[d] > Ubd[d]){
				//Rprintf("Problems with UBD in addPts.\n");
			} else {
				Ubd[d] = UBD[d];
			}
			if(LBD[d] < Lbd[d]){
				//Rprintf("Problems with LBD in addPts.\n");
			} else {
				Lbd[d] = LBD[d];
			}
		}
		lbd == 0;
		for(p=0; p<*nPTC; p++){
			if(PTC[p] < UBD[0] && PTC[p] > LBD[0]){
				lbd++;
			}
		}
		
		// find the lattice
		findLattice(Lbd, Ubd, dim,
					PTC, nPTC, PTCstart, kern, ppd,
					&maxppd, space, L, &nL);
		stop++;
		// if nL == 1, then only one lattice point remaining
		//if(stop > 20){
		//	Rprintf("Warning: problems with number of 'zooms' required. Aborted.\n");
		//}
		if(nL == 1 | stop > 20){
			stillSearching=0;
		}
	}
	
	// add the single lattice point to PT
	for(d=0; d<*dim; d++){
		PT[*dim * *nPT + d] = L[d];
	}
	nPT += 1;
	*newCost = Lmin;
}

void margPT(double *pts, int *npts, int *key,
	    int *dim, int *keyMax, int *maxObs,
	    double *pm,  double *pa,double *pd,
	    double *wts, int *maxBranch,
	    int *ppd, double *space,
	    double *PT,  int *nPT, double *PTcost,
	    double *costEps, double *posEps,
	    int *match, int *euclid,
	    double *lossOrder, int *useKernOnly,
		int *kernProto)
{
	int pt, k, p, d, nPTC, PTCstart[*dim], adding,
		nPTtemp, LMaxElmts, stop;
	double PTC[*npts * *dim], PTtempCost, R,
		PTtemp[*maxObs* *dim], LBD[*dim], UBD[*dim],
		Lbd[*dim], Ubd[*dim];
		// setup LMaxElmts
		LMaxElmts=ppd[0];
	
	for(d=1; d<*dim; d++){
		LMaxElmts = LMaxElmts * ppd[d];
	}
	for(p=0; p<*npts; p++){
		match[p] = -1;
	}
	
	// put together LBD and UBD
	for(d=0; d<*dim; d++){
		LBD[d] = pts[*npts*d];
		UBD[d] = pts[*npts*d];
		for(p=0; p<*npts; p++){
			if(pts[*npts*d + p] < LBD[d]){
				LBD[d] = pts[*npts*d + p];
			}
			if(pts[*npts*d + p] > UBD[d]){
				UBD[d] = pts[*npts*d + p];
			}
		}
		R = UBD[d] - LBD[d];
		UBD[d] += R/100;
		LBD[d] += -R/100;
	}
	
	getPTC(pts, npts, dim,
	       PTC, &nPTC, PTCstart, posEps);
	
	for(p=0; p<*maxObs* *dim; p++){
		PTtemp[p] = -1;
	}
	nPTtemp=0;
	*PTcost=0;
	for(p=0; p<*npts; p++){
		*PTcost += *pa * wts[key[p]-1];
	}
	adding=1;
	stop=0;
	while(adding){
	
		// revise Lbd and Ubd
		for(d=0; d<*dim; d++){
			Lbd[d] = LBD[d];
			Ubd[d] = UBD[d];
		}
		
		// trim points for
		
		// add a prototype point
		addPTPt(pts, npts, key, match,
				keyMax, dim, maxObs,
				PTC, &nPTC, PTCstart,
				pm, pa, pd,
				wts, maxBranch,
				PTtemp, &nPTtemp, costEps,
				ppd, &LMaxElmts, space,
				Lbd, Ubd, &PTtempCost,
				euclid, lossOrder, useKernOnly,
				kernProto);
		// verify we want to add the point
		if(PTtempCost > *PTcost){
			adding = 0;
		} else { // add on to PT
			for(d=0; d<*dim; d++){
				PT[*nPT* *dim+d] = PTtemp[*nPT* *dim+d];
			}
			nPTtemp++; // probably can drop this line
			*nPT += 1;
			*PTcost = PTtempCost;
		}
		adjustMatch(pts, npts, key, match,
					dim, keyMax, maxObs,
					PT, nPT,
					pm, pa, pd, maxBranch,
					wts, costEps, euclid,
					lossOrder, kernProto);
		stop++;
		if(stop > *maxObs && adding != 0){
			adding=0;
			//Rprintf("Manually stopped.\n");
		}
	}
	
	/*
	if(match[0]){
		for(k=0; k<*key; k++){
			// add this option later
		}
	}
	*/
}



/* BELOW IS CODE FOR KERNEL SMOOTHING APPROACH */

void kern(double *pts, int *npts, int *key, int *match,
		  int *keyMax, int *dim, double *pt,
		  double *pm,  double *pa,  double *pd,
		  double *Kern, double *wts,
		  int *euclid, double *lossOrder)
{
	int d, p, k;
	double cost[*keyMax], temp, temp1, hold;
	
	for(k=0; k<*keyMax; k++){
		cost[k] = *pd;
	}
	for(p=0; p<*npts; p++){
		if(match[p] == -1){
			temp = 0;
			temp1 = 0;
			for(d=0; d<*dim; d++){
				if(euclid[d] == 0){
					temp += pm[d] * pow(fabs(pt[d] - pts[d* *npts + p]), *lossOrder);
				} else {
					temp1 += pm[d]*pm[d]*(pt[d]-pts[d* *npts + p])*(pt[d]-pts[d* *npts + p]);
				}
			}
			hold = temp + pow(sqrt(temp1), *lossOrder) - *pa;
			if(hold < cost[key[p]-1]){ // keys start at 1
				cost[key[p]-1] = hold;
			}
		}
	}
	*Kern = 0;
	for(k=0; k<*keyMax; k++){
		*Kern += wts[k]*cost[k];
	}
}

void kernAddPt(double *pts, int *npts, int *key, int *match,
			   int *keyMax, int *dim,
			   double *pm,  double *pa,  double *pd,
			   double *wts,
			   int *euclid, double *lossOrder,
			   double *PT, int *nPT,
			   double *PTC, int *nPTC)
{
	int d, p, k, bestP, nTempPT, matched[*keyMax];
	double bestKern, tempKern, pt[*dim], tempPT[*nPT * *dim];
	double cost[*keyMax], temp, temp1, hold;
	
	bestKern = *npts* *pd;
	bestP = -1;
	// identify best kernel
	for(p=0; p<*nPTC; p++){
		for(d=0; d<*dim; d++){
			pt[d] = PTC[d* *nPTC + p];
		}
		tempKern = 0;
		kern(pts, npts, key, match, keyMax, dim, pt,
			 pm, pa, pd, &tempKern, wts, euclid, lossOrder);
		if(tempKern < bestKern){
			bestKern = tempKern;
			bestP = p;
		}
	}
	
	if(bestKern < 0){ // update prototype
		for(p=0; p<*nPT; p++){
			for(d=0; d<*dim; d++){
				tempPT[*dim*p + d] = PT[*nPT*d + p];
			}
		}
		for(p=0; p<*nPT; p++){
			for(d=0; d<*dim; d++){
				PT[(*nPT+1)*d + p] = tempPT[*dim*p + d];
			}
		}
		*nPT += 1;
		for(d=0; d<*dim; d++){
			PT[*nPT*(d+1) - 1] = PTC[d* *nPTC + bestP];
		}
		
		// update match
		for(k=0; k<*keyMax; k++){
			matched[k] = -1;
		}
		for(k=0; k<*keyMax; k++){
			cost[k] = *pd;
		}
		for(d=0; d<*dim; d++){
			pt[d] = PTC[d* *nPTC + bestP];
		}
		for(p=0; p<*npts; p++){
			if(match[p] == -1){
				temp = 0;
				temp1 = 0;
				for(d=0; d<*dim; d++){
					if(euclid[d] == 0){
						temp += pm[d] * pow(fabs(pt[d] - pts[d* *npts + p]), *lossOrder);
					} else {
						temp1 += pm[d]*pm[d]*(pt[d]-pts[d* *npts + p])*(pt[d]-pts[d* *npts + p]);
					}
				}
				hold = temp + pow(sqrt(temp1), *lossOrder) - *pa;
				if(hold < cost[key[p]-1]){ // keys start at 1
					cost[key[p]-1] = hold;
					matched[key[p]-1] = p;
				}
			}
		}
		hold = 0;
		for(k=0; k<*keyMax; k++){
			if(matched[k] > -1){
				match[matched[k]] = *nPT-1;
			}
			hold += cost[k];
		}
	} else {
		// Could compute distances to see
		// if worth adding point but that
		// was not in original algorithm.
		// This may slow things down.
	}
	
}

void makePTC(double *pts, int *npts, int *dim, double *posEps,
			 double *PTC, int *nPTC, int *ppd, int *ptcSize)
{
	int d, k, p, nptc, keep, ptcStart[*dim+1];
	int inc[*dim], nD[*dim], go;
	double ptc[*ptcSize], temp;
	
	nptc=0;
	for(d=0; d<*dim; d++){
		ptcStart[d] = nptc;
		for(p=0; p<*npts; p++){
			keep = 1;
			for(k=ptcStart[d]; k<nptc; k++){
				temp = fabs(ptc[k] - pts[*npts*d + p]);
				if(temp < posEps[d]){
					keep = 0;
				}
			}
			if(keep == 1){
				ptc[nptc] = pts[d* *npts + p];
				nptc++;
			}
			if(ppd[d] == nptc - ptcStart[d]){
				break;
			}
		}
	}
	ptcStart[*dim] = nptc;
	*nPTC = 1;
	for(d=0; d<*dim; d++){
		nD[d] = ptcStart[d+1]-ptcStart[d];
		*nPTC = *nPTC*nD[d];
	}
	
	for(d=0; d<*dim; d++){
		inc[d]=0;
	}
	p = 0;
	while(inc[0] < nD[0] && p < *nPTC){
		for(d=0; d<*dim; d++){
			PTC[d * *nPTC + p] = ptc[ptcStart[d]+inc[d]];
		}
		p++;
		go=1;
		for(d=*dim-1; d>-1; d--){
			if(go){
				inc[d] += 1;
				if(inc[d] < nD[d]){
					go=0;
				} else if(d > 0){
					inc[d]=0;
				}
			}
		}
	}
}

void kernProto(double *pts, int *npts, int *key, int *match,
			   int *keyMax, int *dim,
			   double *pm,  double *pa,  double *pd,
			   double *wts,
			   int *euclid, double *lossOrder,
			   double *PT, int *nPT,
			   double *posEps, int *ppd,
			   int *maxPTC, double *cost)
{
	int d, k, p, stillAdding, nPTC, ptSize, ptcSize;
	double PTC[*maxPTC * *dim], temp, temp1, hold;
	int matched[*keyMax];
	
	// find PTC & nPTC
	nPTC = *maxPTC;
	ptcSize = 0;
	for(d=0; d<*dim; d++){
		ptcSize += ppd[d];
	}
	makePTC(pts, npts, dim, posEps, PTC, &nPTC, ppd, &ptcSize);
	
	*nPT = 0;
	stillAdding = 1;
	while(stillAdding == 1){
		stillAdding = 0;
		ptSize = *nPT;
		kernAddPt(pts, npts, key, match,
				  keyMax, dim,
				  pm,  pa,  pd,
				  wts,
				  euclid, lossOrder,
				  PT, nPT,
				  PTC, &nPTC);
		if(*nPT > ptSize){
			stillAdding = 1;
		}
	}
	
	for(k=0; k<*keyMax; k++){
		matched[k] = 0;
	}
	*cost = 0;
	for(p=0; p<*npts; p++){
		if(match[p] > -1){
			temp = 0;
			temp1 = 0;
			for(d=0; d<*dim; d++){
				if(euclid[d] == 0){
					temp += pm[d] * pow(fabs(PT[d* *nPT + match[p]] - pts[d* *npts + p]), *lossOrder);
				} else {
					temp1 += pm[d]*pm[d]*pow(fabs(PT[d* *nPT + match[p]]-pts[d* *npts + p]), 2);
				}
			}
			*cost += temp + pow(sqrt(temp1), *lossOrder);
		} else {
			*cost += wts[key[p]] * *pa;
		}
	}
	for(k=0; k<*keyMax; k++){
		if(matched[k] < *nPT){
			*cost += wts[k] * *pd * (*nPT - matched[k]);
		}
	}
}

/////////////////////////////
// VP97 algorithm, supporting
// functions, and method for
// prototypes
/////////////////////////////

void findMin(double x1, double x2, double x3,
			 double *theMin, int *choice)
{
	*theMin = x1;
	*choice = 1;
	if(x2 < x1){
		*theMin = x2;
		*choice = 2;
	}
	if(x3 < *theMin){
		*theMin = x3;
		*choice = 3;
	}
}

void sortThis(double *p, int *np)
{
	int i, j, used[*np], theMin;
	double op[*np];
	
	if(*np > 0){
		theMin = 0;
		for(i=0; i<*np; i++){
			op[i] = p[i];
			if(op[i] < op[theMin]){
				theMin = i;
			}
			used[i] = 0;
		}
		p[0] = op[theMin];
		used[theMin] = 1;
		for(i=1; i<*np; i++){
			theMin = -1;
			for(j=0; j<*np; j++){
				if(used[j] == 0){
					if(theMin == -1){
						theMin = j;
					} else if(op[j] < op[theMin]){
						theMin = j;
					}
				}
			}
			p[i] = op[theMin];
			used[theMin] = 1;
		}
		
	}
}

void vpAlg(double *x, double *y, int *nx, int *ny,
		   double *pm, double *pa, double *pd,
		   int *xtoy, int *ytox, double *cost)
{
	int    i, j, choice, searching;
	double G[*nx+1][*ny+1], g[*nx+1][*ny+1];
	double temp, theMin;
	double t1, t2, t3;
	
	G[0][0] = 0;
	for(i=1; i<*nx+1; i++){
		G[i][0] = G[i-1][0] + *pd;
	}
	for(j=1; j<*ny+1; j++){
		G[0][j] = G[0][j-1] + *pa;
	}
	for(i=1; i<*nx+1; i++){
		for(j=1; j<*ny+1; j++){
			t1 = G[i-1][j]+*pd;
			t2 = G[i][j-1]+*pa;
			t3 = G[i-1][j-1] + *pm*fabs(x[i-1]-y[j-1]);
			findMin(t1, t2, t3, &theMin, &choice);
			G[i][j] = theMin;
			g[i][j] = choice;
		}
	}
	*cost = G[*nx][*ny];
	for(i=0; i<*nx; i++){
		xtoy[i] = -1;
	}
	for(j=0; j<*ny; j++){
		ytox[j] = -1;
	}
	i = *nx;
	j = *ny;
	while(i > 0 && j > 0){
		if(g[i][j] == 3){
			// x[i-1] matches to y[j-1], bounce to i-1, j-1
			xtoy[i-1] = j-1;
			ytox[j-1] = i-1;
			i--;
			j--;
		}
		if(g[i][j] == 2){
			// no matching
			j--;
		}
		if(g[i][j] == 1){
			// no matching
			i--;
		}
	}
}


void vpProto(double *pts, int *npts, int *key, int *keyMax,
			 int *maxObs, double *pm, double *pa, double *pd,
			 double *wts, double *PT, int *nPT, double *ptCost,
			 double *costEps, double *posEps,
			 int *Match, double *lossOrder)
{
	int    i, j, npp;
	int    nPTC, keep, stillLooking=1, best;
	int    npatt, KEY, intTemp, XTOY[*maxObs], YTOX[*maxObs];
	double PTC[*npts], cost, newCost, oldCost, bestCost;
	double patt[*maxObs], hold1, hold2, holdThis;
	double theCost[*keyMax];
	int npt;
	double pt[2* *maxObs];
	
	*ptCost = 0;
	for(KEY=1; KEY<*keyMax+1; KEY++){
		intTemp = 0;
		for(j=0; j<*npts; j++){
			if(key[j] == KEY){
				intTemp++;
			}
		}
		theCost[KEY-1] = *pa * intTemp * wts[KEY-1];
		*ptCost = *ptCost + theCost[KEY-1];
	}
	*nPT = 0;
	PTC[0] = pts[0];
	nPTC = 1;
	for(i=0; i<*npts; i++){
		keep = 1;
		for(j=0; j<nPTC; j++){
			if(fabs(pts[i] - PTC[j]) <= *posEps){
				keep = 0;
			}
		}
		if(keep == 1){
			PTC[nPTC] = pts[i];
			nPTC++;
		}
	}
	
	oldCost = *ptCost;
	
	// look for prototype points
	while(stillLooking == 1 && *nPT < *maxObs){
		stillLooking = 0;
		best = -1;
		bestCost = *npts*(*pa+*pd);
		npt = *nPT+1;
		
		// propose point
		for(npp=0; npp<nPTC; npp++){
			for(i=0; i<*nPT; i++){
				pt[i] = PT[i];
			}
			pt[*nPT] = PTC[npp];
			sortThis(pt, &npt);
			newCost = 0;
			for(KEY=1; KEY<*keyMax+1; KEY++){
				// build pattern
				npatt=0;
				for(i=0; i<*npts; i++){
					if(key[i] == KEY){
						patt[npatt] = pts[i];
						npatt++;
					}
				}
				cost = 0;
				sortThis(patt, &npatt);
				if(npatt == 0){
					cost = *pd * npt;
				} else {
					vpAlg(pt, patt, &npt, &npatt,
						  pm, pa, pd, XTOY, YTOX, &cost);
				}
				newCost = newCost + wts[KEY-1]*cost;
			}
			if(newCost < bestCost - *costEps){
				best = npp;
				bestCost = newCost;
			}
		}
		
		if(bestCost < oldCost - *costEps){
			PT[*nPT] = PTC[best];
			*nPT     = *nPT + 1;
			stillLooking = 1;
			oldCost  = bestCost;
		}
	}
	*ptCost = oldCost;
}


double ppdist093004(int n1, int n2, double *x, double *y, double pen1, double pa) {
	int i,j,k,nx,ny,sumx,sumy;
	double a1,a2;
	double shrt[1000];
	double lng[1000];
	double closest1[1000];
	double closest2[1000];
	int partner1[1000];
	int partner2[1000];
	int keep1[1000];
	int newx[1000];
	int newy[1000];
	int newn1, newn2, stp;
	int i1,i2,i3,i4,i5,i6,i7,i8,i9,keeptot,j1,j2,j3,j4;
	int keepx[1000];
	int keepy[1000];
	double dist,disttemp;
	
	nx = n1;
	ny = n2;
	
	// Rprintf("\n %d %d %f ....",n1,n2,pa);
	// for(i=0;i<n1;i = i+1) Rprintf("%f, ",x[i]);
	// Rprintf(".....");
	// for(i=0;i<n2;i++) Rprintf("%f, ",y[i]);
	// Rprintf("...All done. ......");
	
	
	// // Determine which is the short string and which is longer.
	// I didn't end up using this part.... yet. But might if I do something
	// more efficient later. 	    
	
	// if(nx < ny) {
	//  ns = nx;
	//  nl = ny;
	// for(i=0;i<ns;i++) shrt[i] = x[i];
	//    for(i=0;i<nl;i++) lng[i] = y[i];
	// } 
	// else {
	//	ns = ny;
	//	nl = nx;
	//	for(i=0;i<ns;i++) shrt[i] = y[i];
	//       for(i=0;i<nl;i++) lng[i] = x[i];
	//   }
	
	
	// 1) For each point of x, find the point of y closest to it.
	// And vice versa. This is the "partner" and the distance is "closest".
	
	
	if(nx < 2){
		if(nx == 0){
			dist = pa*ny;
		} else if(ny == 0){
			dist = pa*nx;
		} else {
			a1 = pen1*fabs(x[0] - y[0]);
			for(i=1; i<ny; i++){
				if(a1 > pen1*fabs(x[0] - y[i])){
					a1 = pen1*fabs(x[0] - y[i]);
				}
			}
			if(a1 < 2*pa){
				dist = a1 + (ny-1)*pa;
			} else {
				dist = (ny+1)*pa;
			}
		}
	} else if(ny < 2){
		if(ny == 0){
			dist = pa*nx;
		} else {
			a1 = pen1*fabs(x[0] - y[0]);
			for(i=1; i<nx; i++){
				if(a1 > pen1*fabs(x[i] - y[0])){
					a1 = pen1*fabs(x[i] - y[0]);
				}
			}
			if(a1 < 2*pa){
				dist = a1 + (nx-1)*pa;
			} else {
				dist = (nx+1)*pa;
			}
		}
	} else {
		for(i = 0; i < nx; i++) {
			a1 = fabs(x[i]-y[0]);
			partner1[i] = 0;
			k = 0;
			if(ny>1) {
				for(j=1; k<1; j++) {
					a2 = fabs(x[i]-y[j]);
					if(a2 > a1) k = 2;
					else {		
						a1 = a2;	            
						partner1[i] = j;	        
					}
					if(j == ny-1) k = 2;
				}
			}	    
			closest1[i] = a1;
		}
		
		
		
		for(i = 0; i < ny; i++) {
			a1 = fabs(y[i]-x[0]);
			partner2[i] = 0;
			k = 0;
			if(nx>1) {
				for(j=1; k<1; j++) {	        
					a2 = fabs(y[i]-x[j]);	        
					if(a2 > a1) k = 2;	        
					else {		    
						a1 = a2;	            
						partner2[i] = j;
					}
					if(j == nx-1) k = 2;
				}  
			}
			closest2[i] = a1;
		}
		
		
		//	Rprintf("\n Closest to x are ");
		//	for(i = 0; i < nx; i++) Rprintf(" %f ",closest1[i]);
		//	Rprintf("\n Indices are ");
		//	for(i = 0; i < nx; i++) Rprintf(" %d ",partner1[i]);
		//	Rprintf("\n Closest to y are ");
		//	for(i = 0; i < ny; i++) Rprintf(" %f ",closest2[i]);
		//	Rprintf("\n Indices are ");
		//	for(i = 0; i < ny; i++) Rprintf(" %d ",partner2[i]);
		
		// // 2) If a point in x and y are mutual partners, and if they're closer than 2pa/pen1,
		// then we KNOW they get kept. 
		// So set keep1 to 2 for both of these points.
		// keep1 is a single vector of length nx + ny, representing both x and y.
		// count how many x's you kept this way: the total is keeptot. 
		
		keeptot = 0;
    	for(i = 0; i < nx+ny; i++) keep1[i] = 0;
    	for(i = 0; i < nx; i++) {
    	    if(partner2[partner1[i]] == i) {
    	        if(pen1 * closest1[i] < 2.0 * pa) {
					//    	            Rprintf("\n\n Definitely keep x[%d] and y[%d]",i+1,partner1[i]+1);
					keep1[i] = 2;
					keep1[nx+partner1[i]] = 2;
					keeptot ++;
    	        }
    	    }
    	    if(pen1 * closest1[i] > 2 * pa) {
    	        keep1[i] = 4; 
				//		Rprintf("\n\n Definitely remove x[%d]",i+1);
    	    }
    	}
    	
		// 3) If x (or y) has its partner > 2pa/pen1 away, then we KNOW it gets deleted.
		// So set keep1 to 4 for this point.
		
		for(i = 0; i < ny; i++) {
			if(pen1 * closest2[i] > 2 * pa) {
				keep1[nx+i] = 4;
				//		Rprintf("\n\n Definitely remove y[%d]",i+1);    
			}
		}
		
		
    	
		// 4) Now, count up the remaining points. There'll be newn1 in x and newn2 in y.
		// Let newx be a vector of length newn1, so that newx[0] is the index in x
		// of the first point that might or might not be kept.
		// That is, x[newx[0]] is the first point, x[newx[1]] is the next one, etc.
		
		j = 0;
		for(i = 0; i < nx; i++) {
			if(keep1[i]<1) {
				newx[j] = i;
				j++;
			}   
		}
		newn1 = j;
		
		j = 0;
		for(i = 0; i < ny; i++) {
			if(keep1[nx+i]<1) {
				newy[j] = i;
				j++;
			}
		}
		newn2 = j;
		//	Rprintf("\n newn1 = %d, newn2 = %d, \n",newn1, newn2);
		
		// 5) Now, keepx will be the kept newx's. keepx is a vector of length newn1,
		// starting as all 0's which means none of the points is kept.
		// Store sumx = the number of kept newx's.
		// Same for y.
		// Then go through all possible y's:
		// 	Toggle digit #1. If it goes 0 -> 1, then just increase sumy by 1.
		//	If it goes 1 -> 0, then decrease sumy by 1 and toggle digit #2, etc.
		// 	Stop when you toggle the last digit, i.e. after sumy = newn2.
		//
		// If sumx = sumy, then calculate distance and keep it if it's less than previous. 
		// Then repeat, going through all possible x's the same way.
		
		
		dist = 99999999.9;
		
		for(i = 0; i < newn1; i ++) keepx[i] = 0;
		sumx = 0;
		i4 = 0;
		while(i4 < 1) {
			for(i = 0; i < newn2; i ++) keepy[i] = 0; 
			sumy = 0;
			i1 = 0;
			while(i1 < 1) {  	
				//		Rprintf("\n i = %d, sumx = %d, sumy = %d, ",i,sumx,sumy);	
				//		Rprintf("\n x is ");
				//		for(j4=0;j4<newn1; j4++) Rprintf(" %d ",keepx[j4]);
				//		Rprintf("\n y is ");
				//		for(j4=0;j4<newn2; j4++) Rprintf(" %d ",keepy[j4]);
				//		Rprintf("\n");
				if(sumx == sumy) { 
					disttemp = 0.0;
					// now, calculate the distance:
					// Find the smallest kept x [indexed by j1]
					// Count x[j1] if keep1[j1]=2, or if j1=newx[i8] and keepx[i8]=1
					// and the smallest kept y, [indexed by j3]
					// and then find the distance between them,
					// then repeat. 
					i8 = 0;
					i9 = 0;
					j1 = 0;
					j3 = 0;
					for(i7 = 0; i7 < sumx+keeptot; i7 ++) {
						j2 = 0; // j2 just tells us when to stop looking for next kept x.
						while(j2<1) {
							if(keep1[j1]==2) {
								a1 = x[j1];
								j2 = 2;
							} 
							if(keep1[j1]==0) {
								if(keepx[i8]==1) {
									a1 = x[j1];
									j2 = 2;
								}
								i8++;
							}
							j1 ++;
						}
						
						// Now, find the next y-value.
						j2 = 0; // j2 now tells us when to stop looking for next kept y.
						while(j2<1) {
							if(keep1[nx+j3]==2) {
								a2 = y[j3];
								j2 = 2;
							} 
							if(keep1[nx+j3]==0) {
								if(keepy[i9]==1) {
									a2 = y[j3];
									j2 = 2;
								}
								i9++;
							}
							j3 ++;
						}
						
						//			Rprintf("\n %f %f %f",a1,a2,fabs(a1-a2));
						disttemp += pen1*fabs(a1-a2);
					}
					//		    Rprintf("\nInitial temp distance before penalty is %f\n",disttemp);
					disttemp += pa * (nx+ny-keeptot-keeptot-sumx-sumy);
					//		    Rprintf("X's kept are ");	
					//		    for(i=0;i<newn1; i++) if(keepx[i] == 1) Rprintf(" %f ",x[newx[i]]);	
					//		    for(i=0;i<nx; i++) if(keep1[i] == 2) Rprintf(" %f ",x[i]);	
					//		    Rprintf("\n Y's kept are ");	
					//		    for(i=0;i<newn1; i++) if(keepy[i] == 1) Rprintf(" %f ",y[newy[i]]);	
					//		    for(i=0;i<ny; i++) if(keep1[nx+i] == 2) Rprintf(" %f ",y[i]);	
					//		    Rprintf("\n");
					//		    Rprintf("\nTemp distance = %f,\n\n",disttemp);
					if(disttemp < dist) {
						dist = disttemp;
						//			Rprintf("\n\n distance = %f,\n", dist);
						//			Rprintf("X's kept are ");
						//			for(i=0;i<newn1; i++) if(keepx[i] == 2) Rprintf(" %f ",x[newx[i]]);
						//			for(i=0;i<nx; i++) if(keep1[i] == 2) Rprintf(" %f ",x[i]);
						//			Rprintf("\n Y's kept are ");
						//			for(i=0;i<newn1; i++) if(keepy[i] == 2) Rprintf(" %f ",y[newy[i]]);
						//			for(i=0;i<ny; i++) if(keep1[nx+i] == 2) Rprintf(" %f ",y[i]);
						//			Rprintf("\n");
					}
				}
				if(sumy > newn2-1) i1 = 2;
				i2 = 0;
				i3 = 0;
				while(i2 < 1) {
					if(keepy[i3]==0) {
						keepy[i3]=1;
						sumy += 1;
						i2 = 2;
					}
					else {
						keepy[i3]=0;
						sumy -= 1;
						i3 ++;
						if(i3 == newn2) i2 = 2;
					}
				}
			}
			
			
			
			if(sumx > newn1-1) i4 = 2;
			i5 = 0;
			i6 = 0;
			while(i5 < 1) {
				if(keepx[i6]==0) {
					keepx[i6]=1;
					sumx += 1;
					i5 = 2;
				} 
				else {
					keepx[i6]=0;
					sumx -= 1;
					i6 ++;
					if(i6 == newn1) i5 = 2;
				}
				
			}
			
		}
		// Rprintf("... %f done.",dist);
	}
	return(dist);	
}	

void msuDist(int *n1, int *n2, double *x, double *y, double *pa, double *pm, double *ans2) {
    *ans2 = ppdist093004(*n1,*n2,x,y,*pm,*pa);
}


/* HERE AND BELOW IS NEW (April, 2010) */

void msuProto(double *pts, int *npts, int *key, int *keyMax,
			  int *maxObs, double *pm, double *pa, double *pd,
			  double *wts, double *PT, int *nPT, double *ptCost,
			  double *costEps, double *posEps,
			  int *Match, double *lossOrder)
{
	int    i, j, npp;
	int    nPTC, keep, stillLooking=1, best;
	int    npatt, KEY, intTemp, XTOY[*maxObs], YTOX[*maxObs];
	double PTC[*npts], cost, newCost, oldCost, bestCost;
	double patt[*maxObs], hold1, hold2, holdThis;
	double theCost[*keyMax];
	int npt;
	double pt[2* *maxObs];
	
	*ptCost = 0;
	for(KEY=1; KEY<*keyMax+1; KEY++){
		intTemp = 0;
		for(j=0; j<*npts; j++){
			if(key[j] == KEY){
				intTemp++;
			}
		}
		theCost[KEY-1] = *pa * intTemp * wts[KEY-1];
		*ptCost = *ptCost + theCost[KEY-1];
	}
	*nPT = 0;
	PTC[0] = pts[0];
	nPTC = 1;
	for(i=0; i<*npts; i++){
		keep = 1;
		for(j=0; j<nPTC; j++){
			if(fabs(pts[i] - PTC[j]) <= *posEps){
				keep = 0;
			}
		}
		if(keep == 1){
			PTC[nPTC] = pts[i];
			nPTC++;
		}
	}
	
	oldCost = *ptCost;
	
	// look for prototype points
	while(stillLooking == 1 && *nPT < *maxObs){
		stillLooking = 0;
		best = -1;
		bestCost = *npts*(*pa+*pd);
		npt = *nPT+1;
		
		// propose point
		for(npp=0; npp<nPTC; npp++){
			for(i=0; i<*nPT; i++){
				pt[i] = PT[i];
			}
			pt[*nPT] = PTC[npp];
			sortThis(pt, &npt);
			newCost = 0;
			for(KEY=1; KEY<*keyMax+1; KEY++){
				// build pattern
				npatt=0;
				for(i=0; i<*npts; i++){
					if(key[i] == KEY){
						patt[npatt] = pts[i];
						npatt++;
					}
				}
				cost = 0;
				sortThis(patt, &npatt);
				if(npatt == 0){
					cost = *pd * npt;
				} else {
					msuDist(&npt, &npatt, pt, patt, pa, pm, &cost);
				}
				newCost = newCost + wts[KEY-1]*cost;
			}
			if(newCost < bestCost - *costEps){
				best = npp;
				bestCost = newCost;
			}
		}
		
		if(bestCost < oldCost - *costEps){
			PT[*nPT] = PTC[best];
			*nPT     = *nPT + 1;
			stillLooking = 1;
			oldCost  = bestCost;
		}
	}
	*ptCost = oldCost;
}


