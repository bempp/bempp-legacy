namespace amos
{

int zairy(double *zr, double *zi, int *id, int *kode,
    double *air, double *aii, int *nz, int *ierr);

int zbiry(double *zr, double *zi, int *id, 
	int *kode, double *bir, double *bii, int *ierr);

int zbesh(double *zr, double *zi, double *fnu, 
	int *kode, int *m, int *n, double *cyr, double *cyi,
	int *nz, int *ierr);

int zbesi(double *zr, double *zi, double *fnu, 
	int *kode, int *n, double *cyr, double *cyi, int *nz, int *ierr);

int zbesj(double *zr, double *zi, double *fnu, 
	int *kode, int *n, double *cyr, double *cyi, int *nz, int *ierr);

int zbesk(double *zr, double *zi, double *fnu, 
	int *kode, int *n, double *cyr, double *cyi, int *nz, int *ierr);

int zbesy(double *zr, double *zi, double *fnu, 
	int *kode, int *n, double *cyr, double *cyi, int *nz,
	double *cwrkr, double *cwrki, int *ierr);

} // namespace amos
