#include "Vector.h"

// Matrices

Matrix MatrixMul(Matrix a,Matrix b)
{
        Matrix res;
        for(int i=0;i<16;i++)
        {
                int row=i&3,column=i&12;
                Scalar val=0;
                for(int j=0;j<4;j++) val+=a.a[row+j*4]*b.a[column+j];
                res.a[i]=val;
        }
        return res;
}

Matrix FastMatrixMul(Matrix a,Matrix b)
{
        return(MakeMatrix(a.a[0]*b.a[0]+a.a[4]*b.a[1]+a.a[8]*b.a[2],
                          a.a[0]*b.a[4]+a.a[4]*b.a[5]+a.a[8]*b.a[6],
                          a.a[0]*b.a[8]+a.a[4]*b.a[9]+a.a[8]*b.a[10],
                          a.a[0]*b.a[12]+a.a[4]*b.a[13]+a.a[8]*b.a[14]+a.a[12],

                          a.a[1]*b.a[0]+a.a[5]*b.a[1]+a.a[9]*b.a[2],
                          a.a[1]*b.a[4]+a.a[5]*b.a[5]+a.a[9]*b.a[6],
                          a.a[1]*b.a[8]+a.a[5]*b.a[9]+a.a[9]*b.a[10],
                          a.a[1]*b.a[12]+a.a[5]*b.a[13]+a.a[9]*b.a[14]+a.a[13],

                          a.a[2]*b.a[0]+a.a[6]*b.a[1]+a.a[10]*b.a[2],
                          a.a[2]*b.a[4]+a.a[6]*b.a[5]+a.a[10]*b.a[6],
                          a.a[2]*b.a[8]+a.a[6]*b.a[9]+a.a[10]*b.a[10],
                          a.a[2]*b.a[12]+a.a[6]*b.a[13]+a.a[10]*b.a[14]+a.a[14],

                      0,0,0,1));
}

Matrix FastMatrixInverse(Matrix m)
{
        Matrix res;
        Scalar det=m.a[0]*m.a[5]*m.a[10]-m.a[0]*m.a[6]*m.a[9]+
                  m.a[1]*m.a[6]*m.a[8]-m.a[1]*m.a[4]*m.a[10]+
                          m.a[2]*m.a[4]*m.a[9]-m.a[2]*m.a[5]*m.a[8];
        // singular if det==0
        det=1/det;

        res.a[0]=(m.a[5]*m.a[10]-m.a[6]*m.a[9])*det;
        res.a[4]=-(m.a[4]*m.a[10]-m.a[6]*m.a[8])*det;
        res.a[8]=(m.a[4]*m.a[9]-m.a[5]*m.a[8])*det;

        res.a[1]=-(m.a[1]*m.a[10]-m.a[2]*m.a[9])*det;
        res.a[5]=(m.a[0]*m.a[10]-m.a[2]*m.a[8])*det;
        res.a[9]=-(m.a[0]*m.a[9]-m.a[1]*m.a[8])*det;

        res.a[2]=(m.a[1]*m.a[6]-m.a[2]*m.a[5])*det;
        res.a[6]=-(m.a[0]*m.a[6]-m.a[2]*m.a[4])*det;
        res.a[10]=(m.a[0]*m.a[5]-m.a[1]*m.a[4])*det;

        res.a[12]=-(m.a[12]*res.a[0]+m.a[13]*res.a[4]+m.a[14]*res.a[8]);
        res.a[13]=-(m.a[12]*res.a[1]+m.a[13]*res.a[5]+m.a[14]*res.a[9]);
        res.a[14]=-(m.a[12]*res.a[2]+m.a[13]*res.a[6]+m.a[14]*res.a[10]);

        res.a[3]=res.a[7]=res.a[11]=0;
        res.a[15]=1;

        return(res);
}

