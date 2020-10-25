bool sieveL(double x, double y){


double r_min = 0.0007;
double r_max = 0.001345;


double vert1 = 0.0133;
double vert2 = 0.01331;
double hor1 = 0.00612;
double hor2 = 0.01092;
double hor3 = 0.0048;
double hor4 = 0.00852;
double hor5 = 0.01812;
double hor6 = 0.00306;
double hor7 = 0.0153;
double hor8 = 0.02142;

x = -x;
// Center hole shifted 2mm away from beamline?
//y = y+0.002;
y = -y; //;0.002;

//Center hole
if( x*x + (y)*(y) < r_max*r_max ) return false;


//We will start with the red holes
for(int i = -2; i < 4; i++){

if( ( (x-i*vert1)*(x-i*vert1) + (y)*(y) ) < r_min*r_min && i!=0  ) return false;



//Columns to the right 
if( ( (x-i*vert1)*(x-i*vert1) + (y-hor1)*(y-hor1) )  < r_min*r_min ) return false;
if( ( (x-i*vert1)*(x-i*vert1) + (y-hor2)*(y-hor2) )  < r_min*r_min ) return false;
if(  ( (x-i*vert1)*(x-i*vert1) + (y-hor2-hor3)*(y-hor2-hor3) ) < r_min*r_min ) return false;
if(  ( (x-i*vert1)*(x-i*vert1) + (y-hor2-2*hor3)*(y-hor2-2*hor3) ) < r_min*r_min ) return false;


//these are to the left 
//For RHRS, these are to the right
if ( ( (x-i*vert1)*(x-i*vert1) + (y+hor1)*(y+hor1) ) < r_min*r_min && i != 2  ) return false;
if(  ( (x-i*vert1)*(x-i*vert1) + (y+2*hor1)*(y+2*hor1) ) < r_min*r_min ) return false;
if(  ( (x-i*vert1)*(x-i*vert1) + (y+3*hor1)*(y+3*hor1) ) < r_min*r_min ) return false;

}



if( (x-2*vert1)*(x-2*vert1) + (y+hor1)*(y+hor1) < r_max*r_max ) return false;//other large hole


//Holes asymmetrically place 
if( ( (x+(2*vert1+vert2))*(x+(2*vert1+vert2)) + (y)*(y) ) < r_min*r_min  ) return false;
if( ( (x+(2*vert1+vert2))*(x+(2*vert1+vert2)) + (y-hor1)*(y-hor1) ) < r_min*r_min   ) return false;
if( ( (x+(2*vert1+vert2))*(x+(2*vert1+vert2)) + (y-hor2)*(y-hor2) ) < r_min*r_min   ) return false;
if( ( (x+(2*vert1+vert2))*(x+(2*vert1+vert2)) + (y-hor2-hor3)*(y-hor2-hor3) ) < r_min*r_min   ) return false;
if( ( (x+(2*vert1+vert2))*(x+(2*vert1+vert2)) + (y-hor2-2*hor3)*(y-hor2-2*hor3) ) < r_min*r_min   ) return false;


//LHRS-to the right 
if( ( (x+(2*vert1+vert2))*(x+(2*vert1+vert2)) + (y+hor1)*(y+hor1) ) < r_min*r_min   ) return false;
if( ( (x+(2*vert1+vert2))*(x+(2*vert1+vert2)) + (y+2*hor1)*(y+2*hor1) ) < r_min*r_min   ) return false;
if( ( (x+(2*vert1+vert2))*(x+(2*vert1+vert2)) + (y+3*hor1)*(y+3*hor1) ) < r_min*r_min ) return false;


//Blue--These are asymmetrically placed
//Row above center hole 
if( ( (x-0.00665)*(x-0.00665) + (y+hor6)*(y+hor6) ) < r_min*r_min )return false;
if( ( (x-0.00665)*(x-0.00665) + (y+hor7)*(y+hor7)   )  < r_min*r_min ) return false;
if( ( (x-0.00665)*(x-0.00665) + (y+hor8)*(y+hor8)   )  < r_min*r_min ) return false;
if( ( (x-0.00665)*(x-0.00665) + (y-hor5)*(y-hor5)   )  < r_min*r_min ) return false;
if( ( (x-0.00665)*(x-0.00665) + (y-hor4)*(y-hor4)   )  < r_min*r_min ) return false;

//Row below center
if( ( (x+0.00665)*(x+0.00665) + (y+hor6)*(y+hor6) ) < r_min*r_min )return false;
if( ( (x+0.00665)*(x+0.00665) + (y+hor7)*(y+hor7)   )  < r_min*r_min ) return false;
if( ( (x+0.00665)*(x+0.00665) + (y+hor8)*(y+hor8)   )  < r_min*r_min ) return false;
if( ( (x+0.00665)*(x+0.00665) + (y-hor4)*(y-hor4)   )  < r_min*r_min ) return false;
if( ( (x+0.00665)*(x+0.00665) + (y-hor5)*(y-hor5)   )  < r_min*r_min ) return false;


if( ( (x+0.00665+vert1)*(x+0.00665+vert1) + (y-hor4)*(y-hor4)   ) < r_max*r_max ) return false;//Larger blue hole



for(int i = 1; i < 3; i++) {

//Two rows from top of sieve
if( ( (x-0.00665-i*vert1)*(x-0.00665-i*vert1) + (y+hor6)*(y+hor6) ) < r_min*r_min ) return false;
if( ( (x-0.00665-i*vert1)*(x-0.00665-i*vert1) +  (y+hor7)*(y+hor7)   )  < r_min*r_min ) return false;
if( ( (x-0.00665-i*vert1)*(x-0.00665-i*vert1) +  (y+hor8)*(y+hor8)   )  < r_min*r_min ) return false;
if( ( (x-0.00665-i*vert1)*(x-0.00665-i*vert1) +  (y-hor4)*(y-hor4)   )  < r_min*r_min ) return false;
if( ( (x-0.00665-i*vert1)*(x-0.00665-i*vert1) +  (y-hor5)*(y-hor5)   )  < r_min*r_min ) return false;


if(i == 1){//Penultimate blue row--center hole is separate
if( ( (x+0.00665+i*vert1)*(x+0.00665+i*vert1) + (y+hor6)*(y+hor6) ) < r_min*r_min ) return false;
if( ( (x+0.00665+i*vert1)*(x+0.00665+i*vert1) +  (y+hor7)*(y+hor7)   )  < r_min*r_min ) return false;
if( ( (x+0.00665+i*vert1)*(x+0.00665+i*vert1) +  (y+hor8)*(y+hor8)   )  < r_min*r_min ) return false;
if( ( (x+0.00665+i*vert1)*(x+0.00665+i*vert1) +  (y-hor5)*(y-hor5)   )  < r_min*r_min ) return false;
}

}

if( ( (x+0.03395)*(x+0.03395) + (y+hor6)*(y+hor6) ) < r_min*r_min ) return false;

if( ( (x+0.03395)*(x+0.03395) +  (y+hor7)*(y+hor7)   )  < r_min*r_min ) return false;
if( ( (x+0.03395)*(x+0.03395) +  (y+hor8)*(y+hor8)   )  < r_min*r_min ) return false;
if( ( (x+0.03395)*(x+0.03395) +  (y-hor5)*(y-hor5)   )  < r_min*r_min ) return false;
if( ( (x+0.03395)*(x+0.03395) + (y-hor4)*(y-hor4)   ) < r_min*r_min ) return false;



return true;



}

