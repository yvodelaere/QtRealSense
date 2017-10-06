#ifndef LDA_INCLUDED
#define LDA_INCLUDED

#include <cv.h>
namespace utw3dface {

class LDA
{
public:
	CvMat *meanZ;
	CvMat *T;
	CvMat *Tinv;
	CvMat *Rho1;
	CvMat *Ddiag;
	CvMat *T1;
	CvMat *S1diag;

	CvMat *Rho1_T;

	int get_means_and_variations(CvMat *Z,CvMat *iZ,CvMat *&meanZ,CvMat *&meanClasses,CvMat *&Zzeromean,CvMat *&Zvariations);

	// robust version
	int get_means_and_variationsR(CvMat *Z,CvMat *iZ,CvMat *&meanZ,CvMat *&meanClasses,CvMat *&Zzeromean,CvMat *&Zvariations);

	void cleanup();
public:
	LDA();
	~LDA();
	int train(CvMat *Z,CvMat *iZ,int nPCA,int nLDA);
	int retrainLDA(CvMat *Z,CvMat *iZ,int nPCA,int nLDA);
	// robust version
	int trainR(CvMat *Z,CvMat *iZ,int nPCA,int nLDA);
	int load(const char *basename);
	int save(const char *basename);
	int load(FILE *f);
	int save(FILE *f);
	CvMat *enroll(CvMat *v);
	CvMat *features(CvMat *v);
	void enroll(double *v,double *v_enrolled);
	void features(double *v,double *v_features);
	double likelihoodratio(CvMat *reference,CvMat *test);
	double likelihoodratio(double *reference,double *test);
	int featurelength() { return T->rows; }
	
	double DIFS(CvMat *v);
	double DFFS(CvMat *v);
};
}

#endif // LDA_INCLUDED
