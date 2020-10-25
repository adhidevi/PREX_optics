bool UpPlaneL(double xup1, double yup1, double xup2, double yup2){

double x1 = 0.0488;
double x2 = 0.0707;
double y1 = -0.0258;
double y2 = 0.0255;


if( xup1 < x1 ) return false;

if( xup2 > x2 ) return false;

if( yup1 < y1 ) return false;

if( yup1 > y2) return false; 


return true;

}
