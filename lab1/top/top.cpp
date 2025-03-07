#include <set>
#include <gmsh.h>
#include <cmath>
#include <iostream>
int main(int argc, char **argv){
	gmsh::initialize();
	gmsh::model::add("t2");
	double pi = 3.14159265;
	double r = 1;
	double R = 5;
	double x1 = R;
	double z1 = 0;
	for (int j = 0; j < 10; j++){
		x1 = R*cos(2*pi*j/10);
		z1 = R*sin(2*pi*j/10);
	for (int i = 0; i < 10; i++){
		gmsh::model::geo::addPoint(x1+(x1/pow((x1*x1+z1*z1),0.5))*r*cos(2*pi*i/10),r*sin(2*pi*i/10),z1+(z1/pow((x1*x1+z1*z1),0.5))*r*cos(2*pi*i/10),0.5,10*j+i);
		gmsh::model::geo::addPoint(x1+(x1/pow((x1*x1+z1*z1),0.5))*3*r/4*cos(2*pi*i/10),3*r/4*sin(2*pi*i/10),z1+(z1/pow((x1*x1+z1*z1),0.5))*3*r/4*cos(2*pi*i/10),0.5,10*j+i+10000);
		if (i>0){
		gmsh::model::geo::addLine(10*j+i-1,10*j+i,10*j+i-1);
		gmsh::model::geo::addLine(10*j+i-1+10000,10*j+i+10000,10*j+i-1+10000);
		}
	} gmsh::model::geo::addLine(10*j+9,10*j,10*j+9);
	gmsh::model::geo::addLine(10*j+9+10000,10*j+10000,10*j+9+10000);
//	gmsh::model::geo::addCurveLoop({10*j,10*j+1,10*j+2,10*j+3,10*j+4,10*j+5,10*j+6,10*j+7,10*j+8,10*j+9},j);
	}
	for (int j = 0; j < 10; j++){
	for (int i = 0; i < 10; i++){
			gmsh::model::geo::addLine(10*j+i,10*((j+1)%10)+i,100+10*j+i);
			gmsh::model::geo::addLine(10*j+i,10*((j+1)%10)+(i+1)%10,200+10*j+i);
			gmsh::model::geo::addLine(10*j+i+10000,10*((j+1)%10)+i+10000,100+10*j+i+10000);
			gmsh::model::geo::addLine(10*j+i+10000,10*((j+1)%10)+(i+1)%10+10000,200+10*j+i+10000);
	}
	}
	for (int j = 0; j < 10; j++){
	for (int i = 0; i < 10; i++){
		gmsh::model::geo::addCurveLoop({100+10*j+i,10*((j+1)%10)+i,-(200+10*j+i)},10*j+i+1);
		gmsh::model::geo::addPlaneSurface({10*j+i+1},10*j+i+1);
		gmsh::model::geo::addCurveLoop({100+10*j+i+10000,10*((j+1)%10)+i+10000,-(200+10*j+i+10000)},10*j+i+1+10000);
		gmsh::model::geo::addPlaneSurface({10*j+i+1+10000},10*j+i+1+10000);
	}
	}
        for (int j = 0; j < 10; j++){
	for (int i = 0; i < 10; i++){
	gmsh::model::geo::addCurveLoop({-(200+10*j+i),(100+10*j+(i+1)%10),(10*j+i)},100+10*j+i+1);
	gmsh::model::geo::addPlaneSurface({100+10*j+i+1},100+10*j+i+1);
	gmsh::model::geo::addCurveLoop({-(200+10*j+i+10000),(100+10*j+(i+1)%10+10000),(10*j+i+10000)},100+10*j+i+1+10000);
	gmsh::model::geo::addPlaneSurface({100+10*j+i+1+10000},100+10*j+i+1+10000);
	}
	}
	std::vector<int> sur1;
	std::vector<int> sur2;
	for (int i = 0; i < 200; i++){
		sur1.push_back(i+1);
		sur2.push_back(i+10001);
	}
	gmsh::model::geo::addSurfaceLoop(sur1,1);
	gmsh::model::geo::addSurfaceLoop(sur2,2);
	gmsh::model::geo::addVolume({1,-2});
	gmsh::model::geo::synchronize();
	gmsh::model::mesh::generate(3);
	gmsh::write("t2.msh");
	std::set<std::string> args(argv, argv + argc);
	if(!args.count("-nopopup")) gmsh::fltk::run();

	gmsh::finalize();

	return 0;
}
