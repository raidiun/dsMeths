#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>

#define PI 3.141592

struct Vector {
	int size;
	double* data;
	};

struct Vector makeVectorOfSize(int size) {
	struct Vector outVector;
	outVector.size = size;
	outVector.data = malloc(sizeof(double)*size);
	return(outVector);
	}

struct Vector makeVector(int size, ...) {
	va_list args;
	int idx;
	
	struct Vector outVector=makeVectorOfSize(size);
	
	va_start(args,size);
	for(idx=0;idx<size;idx++) {
		outVector.data[idx] = va_arg(args,double);
		}
	va_end(args);
	
	return(outVector);
	}

struct Vector makeVectorWithArray(int size, double* array) {
	int idx;
	
	struct Vector outVector=makeVectorOfSize(size);
	
	for(idx=0;idx<size;idx++) {
		outVector.data[idx] = array[idx];
		}
	
	return(outVector);
	}

void vectorPrint(char name[8],struct Vector inVector) {
	int idx=0;
	
	printf("%8s:[",name);
	
	for(idx=1;idx<inVector.size;idx++) {
		printf("%lf,",inVector.data[idx-1]);
		}
	printf("%lf]\n",inVector.data[idx-1]);
	}

double vectorDotProduct(struct Vector a, struct Vector b) {
	int idx;
	double result=0;
	
	if(a.size == b.size) {
		for(idx=0;idx<a.size;idx++) {
			result += a.data[idx] * b.data[idx];
			}
		}
	
	return(result);
	}

double vectorMagnitude(struct Vector inVec) {
	int idx;
	double squaresSum=0;
	for(idx=0;idx<inVec.size;idx++) {
		squaresSum += pow(inVec.data[idx],2);
		}
	return(sqrt(squaresSum));
	}

struct Vector vectorNormalise(struct Vector inVec) {
	int idx;
	double vMag = vectorMagnitude(inVec);
	
	for(idx=0;idx<inVec.size;idx++) {
		inVec.data[idx] /= vMag;
		}
	
	return(inVec);
	}

struct Vector vectorvectorAdd(struct Vector a, struct Vector b) {
	if(a.size != b.size) {
		printf("Vector Add Error: Different Sizes");
		return(makeVectorOfSize(0));
		}
	else {
		struct Vector outVector=makeVectorOfSize(a.size);
		int idx;
		for(idx=0;idx<a.size;idx++) {
			outVector.data[idx] = a.data[idx] + b.data[idx];
			}
		return(outVector);
		}
	}

struct Vector vectorvectorSubtract(struct Vector a, struct Vector b) {//a-b
	struct Vector tVect=makeVectorOfSize(b.size);
	int idx;
	for(idx=0;idx<b.size;idx++) {
		tVect.data[idx] = -1.0*b.data[idx];
		}
	
	return(vectorvectorAdd(a,tVect));
	}

struct Vector vectorScalarMultiply(double scalar, struct Vector inVector) {
	int idx;
	struct Vector outVector=makeVectorOfSize(inVector.size);
	for(idx=0;idx<inVector.size;idx++) {
		outVector.data[idx] = scalar*inVector.data[idx];
		}
	return(outVector);
	}

struct Matrix {
	int size;
	double** data;
	};

struct Matrix makeMatrixOfSize(int size) {
	struct Matrix outMatrix;
	outMatrix.size = size;
	
	outMatrix.data = malloc(sizeof(double*)*outMatrix.size);
	for(;size>0;size--) {
		outMatrix.data[size-1] = malloc(sizeof(double)*outMatrix.size);
		}
	return(outMatrix);
	}

struct Matrix makeMatrix(int size, ...) {
	va_list args;
	int idx=0, sizeSqd=pow(size,2);
	
	struct Matrix outMatrix=makeMatrixOfSize(size);
	
	va_start(args,size);
	for(idx=0;idx<sizeSqd;idx++) {
		outMatrix.data[idx/size][idx%size] = va_arg(args,double);
		}
	va_end(args);
	
	return(outMatrix);
	}

struct Matrix makeMatrixWithArray(int size, double** array) {
	int idx=0, sizeSqd=pow(size,2);
	
	struct Matrix outMatrix=makeMatrixOfSize(size);
	
	for(idx=0;idx<sizeSqd;idx++) {
		outMatrix.data[idx/size][idx%size] = array[idx/size][idx%size];
	}

	return(outMatrix);
}

void matrixPrint(char name[8],struct Matrix inMatrix) {
	int ridx,cidx;
	
	printf("%8s:[",name);
	
	for(ridx=0;ridx<inMatrix.size;ridx++) {
		
		if(ridx!=0) {//If not first row
			printf("          ");
			}
		
		printf("[%f,",inMatrix.data[ridx][0]);
		for(cidx=2;cidx<inMatrix.size;cidx++) {
			printf("%f,",inMatrix.data[ridx][cidx-1]);
			}
		
		if(ridx==(inMatrix.size-1)) {//If last row
			printf("%f]]\n",inMatrix.data[ridx][cidx-1]);
			}
		else {
			printf("%f],\n",inMatrix.data[ridx][cidx-1]);
			}
		}
	}

struct Matrix matrixTranspose(struct Matrix inMatrix) {
	int ridx,cidx,size=inMatrix.size;
	struct Matrix outMatrix=makeMatrixOfSize(size);
	
	for(ridx=0;ridx<size;ridx++) {
		for(cidx=0;cidx<size;cidx++) {
			outMatrix.data[cidx][ridx] = inMatrix.data[ridx][cidx];
			}
		}
	
	return(outMatrix);
	}

struct Matrix matrixSubmatrix(struct Matrix inMatrix,int row,int col) {
	struct Matrix outMatrix=makeMatrixOfSize(inMatrix.size-1);
	
	int irdx=0,icdx=0;
	int ridx,cidx;
	for(ridx=0;ridx<inMatrix.size;ridx++) {
		if(ridx != row) {
			icdx=0;
			for(cidx=0;cidx<inMatrix.size;cidx++) {
				if(cidx != col) {
					outMatrix.data[irdx][icdx] = inMatrix.data[ridx][cidx];
					icdx++;
					}
				}
			irdx++;
			}
		}
	
	return(outMatrix);
	}

double matrixDeterminant(struct Matrix inMatrix) {
	if(inMatrix.size == 2) {
		return((inMatrix.data[0][0]*inMatrix.data[1][1])-(inMatrix.data[0][1]*inMatrix.data[1][0]));
		}
	if(inMatrix.size == 1) {
		return(inMatrix.data[0][0]);
		}
	else {
		int cidx;
		double cofactorness,sum=0;
		for(cidx=0;cidx<inMatrix.size;cidx++) {
			if(cidx%2) {cofactorness = -1.0;} else {cofactorness = 1.0;}
			sum += inMatrix.data[0][cidx]*cofactorness*matrixDeterminant(matrixSubmatrix(inMatrix,0,cidx));
			}
		return(sum);
		}
	}

struct Matrix matrixCofactors(struct Matrix inMatrix) {
	int ridx,cidx,size=inMatrix.size;
	struct Matrix outMatrix=makeMatrixOfSize(size);
	
	double cofactorness;
	for(ridx=0;ridx<size;ridx++) {
		for(cidx=0;cidx<size;cidx++) {
			if((cidx+ridx)%2) {cofactorness = -1.0;} else {cofactorness = 1.0;}
			outMatrix.data[ridx][cidx] = cofactorness*matrixDeterminant(matrixSubmatrix(inMatrix,ridx,cidx));
			}
		}
	return(outMatrix);
	}

struct Matrix matrixInvert(struct Matrix inMatrix) {
	struct Matrix outMatrix;
	
	outMatrix = matrixTranspose(matrixCofactors(inMatrix));
	double det = matrixDeterminant(inMatrix);
	
	int ridx,cidx,size=inMatrix.size;
	for(ridx=0;ridx<size;ridx++) {
		for(cidx=0;cidx<size;cidx++) {
			outMatrix.data[ridx][cidx] = (outMatrix.data[ridx][cidx])/det;
			}
		}
	
	return(outMatrix);
	}

struct Vector matrixvectorMultiply(struct Matrix inMatrix, struct Vector inVector) {
	if(inVector.size != inMatrix.size) {
		printf("Multiply Error: Incompatible Sizes");
		struct Vector outVector;
		outVector.size = 0;
		outVector.data = NULL;
		return(outVector);
		}
	
	else {
		int ridx,cidx,size=inVector.size;
		double vectVals[inVector.size];
		double rowSum;
		
		for(ridx=0;ridx<size;ridx++) {
			rowSum=0;
			for(cidx=0;cidx<size;cidx++) {
				rowSum += (inMatrix.data[ridx][cidx]*inVector.data[cidx]);
				}
			vectVals[ridx] = rowSum;
			}
		
		return(makeVectorWithArray(size,vectVals));
		}
	}

struct Matrix matrixmatrixMultiply(struct Matrix a, struct Matrix b) {
	struct Matrix outMatrix;
	int ridx,cidx,eidx,size=a.size;
	double sum=0;
	
	if(a.size != b.size) {
		printf("Multiply Error: Matricies not of same size");
		outMatrix.size = 0;
		return(outMatrix);
		}
	
	outMatrix=makeMatrixOfSize(size);
	
	for(ridx=0;ridx<size;ridx++) {
		for(cidx=0;cidx<size;cidx++) {
			sum=0.0;
			for(eidx=0;eidx<size;eidx++) {
				sum += a.data[ridx][eidx] * b.data[eidx][cidx];
				}
			outMatrix.data[ridx][cidx] = sum;
			}
		}
	
	return(outMatrix);
	}

struct ControlPoint {
	struct Vector position;
	struct Vector normal;
	};

struct Source {
	struct Vector position;
	double strength;
	};

struct Params {
	int count;
	struct ControlPoint* points;
	struct Source* sources;
	struct Vector freestream;
	};

double influenceOfSourceAtControlPoint(struct Source source,struct ControlPoint ctrl) {
	struct Vector deltaPos=vectorvectorSubtract(ctrl.position,source.position);
	deltaPos = vectorScalarMultiply(1.0/(2.0*PI*pow(vectorMagnitude(deltaPos),2)),deltaPos);
	return(vectorDotProduct(deltaPos,ctrl.normal));
	}

struct Matrix constructInfluenceMatrix(int count, struct Source sources[], struct ControlPoint points[]) {
	int ridx,cidx;
	struct Matrix outMatrix=makeMatrixOfSize(count);
	for(ridx=0;ridx<count;ridx++) {
		for(cidx=0;cidx<count;cidx++) {
			outMatrix.data[ridx][cidx] = influenceOfSourceAtControlPoint(sources[cidx],points[ridx]);
			}
		}
	return(outMatrix);
	}

struct Vector constructFreestreamVector(int count, struct Vector freestream, struct ControlPoint points[]) {
	int idx;
	struct Vector outVector=makeVectorOfSize(count);
	for(idx=0;idx<count;idx++) {
		outVector.data[idx] = -1.0*vectorDotProduct(points[idx].normal,freestream);
		}
	return(outVector);
	}

struct Params parseProblem(int argc, char** argv) {
	int idx,pidx=0,cFlag=0;
	struct Params outParams;
	for(idx=1;idx<argc;) {
		switch(argv[idx][1]) {
			case('f'): {
				struct Vector freestreamVector = makeVectorOfSize(3);
				sscanf(argv[idx+1],"{%lf,%lf,%lf}",&freestreamVector.data[0],&freestreamVector.data[1],&freestreamVector.data[2]);
				outParams.freestream = freestreamVector;
				idx+=2;
				break;
				}
			case('c'): {
				sscanf(argv[idx+1],"%i",&outParams.count);
				outParams.points = malloc(sizeof(struct ControlPoint)*outParams.count);
				outParams.sources = malloc(sizeof(struct Source)*outParams.count);
				cFlag=1;
				idx+=2;
				break;
				}
			case('p'): {
				if(cFlag==0) {
					printf("Problem size not set. Set size before defining control points");
					break;
					}
				if(pidx == outParams.count) {
					printf("Problem full. Increase size\n");
					break;
					}
				
				struct Source source;
				struct ControlPoint point;
				struct Vector pointPos=makeVectorOfSize(3),pointNorm=makeVectorOfSize(3),sourcePos=makeVectorOfSize(3);
				int eidx=1;
				for(;eidx<4;eidx++) {
					switch(argv[idx+eidx][0]) {
						case('s'): {
							sscanf(argv[idx+eidx],"s{%lf,%lf,%lf}",&sourcePos.data[0],&sourcePos.data[1],&sourcePos.data[2]);
							}
						case('n'): {
							sscanf(argv[idx+eidx],"n{%lf,%lf,%lf}",&pointNorm.data[0],&pointNorm.data[1],&pointNorm.data[2]);
							pointNorm = vectorNormalise(pointNorm);
							}
						case('p'): {
							sscanf(argv[idx+eidx],"n{%lf,%lf,%lf}",&pointPos.data[0],&pointPos.data[1],&pointPos.data[2]);
							}
						}
					}
				source.position = sourcePos;
				point.position = pointPos;
				point.normal = pointNorm;
				
				outParams.points[pidx] = point;
				outParams.sources[pidx] = source;
				
				pidx++;
				idx+=4;
				break;
				}
			}
		}
	if(outParams.count != pidx) {
		printf("WARNING: Specified count less than those provided.\n");
		}
	return(outParams);
	}

void printParams(struct Params problem) {
	printf("Problem Description:\nSource Count: %i\n",problem.count);
	vectorPrint("freestr",problem.freestream);
	}

int main(int argc, char** argv) {
	struct Params problem=parseProblem(argc,argv);
	printParams(problem);
	struct Matrix influenceMatrix = constructInfluenceMatrix(problem.count,problem.sources,problem.points);
	struct Vector freestreamVector = constructFreestreamVector(problem.count,problem.freestream,problem.points);
	
	matrixPrint("inflMat",influenceMatrix);
	vectorPrint("fstrVec",freestreamVector);
	
	printf("det(inflMat): %lf\n",matrixDeterminant(influenceMatrix));
	
	influenceMatrix = matrixInvert(influenceMatrix);
	matrixPrint("inflInv",influenceMatrix);
	
	struct Vector sourceStrengths = matrixvectorMultiply(influenceMatrix,freestreamVector);
	
	vectorPrint("srcVec",sourceStrengths);
	/* Steps:
	 - Define Control Points and Normals (Geometry) Y
	 - Define Freestream Y
	 - Define sources (or do this programatically) Y
	 - Create influence matrix
	 - Create freeflow vector
	 - Solve for source strengths
	 - Allow query of values at specified position vector
	 
	 inv(inflMat)*freeVec = sourceStrengths
	 
	 */
	}
