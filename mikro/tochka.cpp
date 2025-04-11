#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <algorithm>

#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkGlyph3D.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkPoints.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkAppendPolyData.h>

using namespace std;
void writeSpheresAsPoints(vtkSmartPointer<vtkPolyData> polyData, const std::string& fileName) {
    // Создаем источник сфер
    vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
    sphereSource->SetRadius(0.01);  // Установливаем радиус сферы
    sphereSource->SetPhiResolution(20);  // Разрешение по углу phi
    sphereSource->SetThetaResolution(20); // Разрешение по углу theta
    sphereSource->Update();

    // Используем Glyph3D для замены точек на сферы
    vtkSmartPointer<vtkGlyph3D> glyph = vtkSmartPointer<vtkGlyph3D>::New();
    glyph->SetSourceConnection(sphereSource->GetOutputPort());
    glyph->SetInputData(polyData);
    glyph->ScalingOff();  
    glyph->Update();
	

    vtkSmartPointer<vtkAppendPolyData> sum_data = vtkSmartPointer<vtkAppendPolyData>::New();
    sum_data->AddInputData(polyData);
    sum_data->AddInputData(glyph->GetOutput());
    sum_data->Update();
    // Создаем writer для записи файла .vtp
    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(fileName.c_str());
    writer->SetInputData(sum_data->GetOutput());
    writer->Write();
}

// Класс расчётной точки
class CalcNode {
public:
    // Координаты
    double x;
    double y;
    double z;
    // Некая величина, в попугаях
    double smth;
    // Скорость
    double vx;
    double vy;
    double vz;

    // Конструктор по умолчанию
    CalcNode() : x(0.0), y(0.0), z(0.0), smth(0.0), vx(0.0), vy(0.0), vz(0.0) {}

    // Конструктор с указанием всех параметров
    CalcNode(double x, double y, double z, double smth, double vx, double vy, double vz) 
        : x(x), y(y), z(z), smth(smth), vx(vx), vy(vy), vz(vz) {}

    // Метод отвечает за перемещение точки
    // Движемся время tau из текущего положения с текущей скоростью
    void move(double tau) {
        x += vx * tau;
        y += vy * tau;
        z += vz * tau;
    }
};
// функция расстояния между точками
double distant(CalcNode _1, CalcNode _2) {
    return sqrt(pow(_1.x - _2.x, 2) + pow(_1.y - _2.y, 2) + pow(_1.z - _2.z, 2));
}

// Класс расчётной сетки
class CalcMesh {
protected:
    // 3D-сетка из расчётных точек
    vector<CalcNode> points;

public:
    // Конструктор сетки size x size точек с шагом h по пространству
    CalcMesh(int N) {
        const unsigned int seed = 228;
        std::mt19937 gen(seed);
        std::uniform_real_distribution<double> distrib(-1.0, 1.0); // Лучше задать реальные пределы
	//главная птичка
	double x = 0;
	double y = 0;
	double z = 0;
	double smth = 1.0;
	double vx = 0;
	double vy = 0;
	double vz = 0;						       
	points.push_back(CalcNode(x, y, z, smth, vx, vy, vz));
	// остальные птички
        for (int i = 1; i < N; i++) {
            double x = distrib(gen);
            double y = distrib(gen);
            double z = distrib(gen);
            double smth = 0.0; // Начальное значение
            double vx = distrib(gen)*0.3;
            double vy = distrib(gen)*0.3;
            double vz = distrib(gen)*0.3;

            points.push_back(CalcNode(x, y, z, smth, vx, vy, vz));
        }
    }

    // Метод отвечает за выполнение для всей сетки шага по времени величиной tau
    void doTimeStep(double tau, double T) {
        double ax;
        double az;
        double ay;
	double abs_v;
	double max_v = 1;
	//тут задаем как летит главная птичка
	if (T < 10){
		points[0].move(tau);
		points[0].vx=0.0;
		points[0].vy=-0.0;
		points[0].vz=0.0;
	}
	if (T >= 10 and T < 17.5){
	points[0].move(tau);
	points[0].vx=0.15;
	points[0].vy=0;
	points[0].vz=0;
	}
	if (T >= 17.5 and T < 28){
	points[0].move(tau);
	points[0].vx=-0.15*(points[0].y-1);
	points[0].vy=0.15*(points[0].x-1);
	points[0].vz=0;
	}
	if (T >=28){
                points[0].move(tau);
		points[0].vx=0.0;
		points[0].vy=0.15;
		points[0].vz=0.0;

	}
	// а дальше остальных
        for(int i = 1; i < points.size(); i++) {
            if (points[i].smth > -0.1){
	    points[i].move(tau);
	    if ((points[i].x>2) and(points[i].y<0)){
		points[i].smth = -0.5;
		if(points[i].x-0.06>2){
			points[i].y = 0.001;}
		if(points[i].y+0.06<0){
			points[i].x = 2 - 0.001;}
		continue;
	    }
            int sosedi[5] = {-1, -1, -1, -1, -1};
            double rasst_sosedi[5] = {1000, 1000, 1000, 1000, 1000};
            for (int k = 0; k < 5; k++){
                for(int j = 0; j < points.size(); j++){
                    if ((j == i) or (find(begin(sosedi), end(sosedi), j)!=end(sosedi))){
                        continue;
                    }
                    if (distant(points[i], points[j]) < rasst_sosedi[k]){
                        rasst_sosedi[k] = distant(points[i], points[j]);
                        sosedi[k] = j;
                    }
                } 
            }
            ax = 0;
            ay = 0;
            az = 0;
            for (int k = 0; k < 5; k++){
                ax += 0.0*points[k].vx + 0.2*(points[k].x-points[i].x)*(points[k].smth+0.25)*4;
                ay += 0.0*points[k].vy + 0.2*(points[k].y-points[i].y)*(points[k].smth+0.25)*4;
                az += 0.0*points[k].vz + 0.3*(points[k].z-points[i].z)*(points[k].smth+0.25)*4;
            }
            points[i].vx += ax * tau;
            points[i].vy += ay * tau;
            points[i].vz += az * tau;
	    abs_v = sqrt(pow(points[i].vx, 2) + pow(points[i].vy, 2) + pow(points[i].vz, 2));
	    if (abs_v > max_v){
		 points[i].vx *= max_v/abs_v;
		 points[i].vy *= max_v/abs_v;
		 points[i].vz *= max_v/abs_v;
	    }
            }
	}
    }

    // Метод отвечает за запись текущего состояния сетки в снапшот в формате VTK
    void snapshot(unsigned int snap_number) {
        // Сетка в терминах VTK
        vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();

        // Точки сетки в терминах VTK
        vtkSmartPointer<vtkPoints> dumpPoints = vtkSmartPointer<vtkPoints>::New();

        // Скалярное поле на точках сетки
        auto smth = vtkSmartPointer<vtkDoubleArray>::New();
        smth->SetName("smth");

        // Векторное поле на точках сетки
        auto vel = vtkSmartPointer<vtkDoubleArray>::New();
        vel->SetName("velocity");
        vel->SetNumberOfComponents(3);
	
        // Обходим все точки нашей расчётной сетки
        int number = (int)points.size();
        for (unsigned int i = 0; i < number; i++) {
            // Вставляем новую точку в сетку VTK-снапшота
            dumpPoints->InsertNextPoint(points[i].x, points[i].y, points[i].z);

            // Добавляем значение векторного поля в этой точке
            double _vel[3] = {points[i].vx, points[i].vy, points[i].vz};
            vel->InsertNextTuple(_vel);

            // И значение скалярного поля тоже
            smth->InsertNextValue(points[i].smth);
        }
        // Создаем препятствие
	vtkSmartPointer<vtkCellArray> cubePolys = vtkSmartPointer<vtkCellArray>::New();
	
	dumpPoints->InsertNextPoint(2.0, 0.0, -2.0);
	dumpPoints->InsertNextPoint(2.0, 0.0, 2.0);
	dumpPoints->InsertNextPoint(4.0, 0.0, -2.0);
	dumpPoints->InsertNextPoint(4.0, 0.0, 2.0);
	dumpPoints->InsertNextPoint(2.0, -2.0, -2.0);
	dumpPoints->InsertNextPoint(2.0, -2.0, 2.0);
	for (int i = 0; i<6; i++){
		double _vel[3] = {0.0, 0.0, 0.0};
		vel->InsertNextTuple(_vel);
		smth->InsertNextValue(-1.0);
	}
	cubePolys->InsertNextCell({0+number, 1+number, 3+number, 2+number});
	cubePolys->InsertNextCell({0+number, 1+number, 5+number, 4+number});
        
	//добавляем точки и полигоны
	polyData->SetPoints(dumpPoints);
	polyData->SetPolys(cubePolys);
        // Присоединяем векторное и скалярное поля к точкам
        polyData->GetPointData()->AddArray(vel);
        polyData->GetPointData()->AddArray(smth);

        // Создаём снапшот в файле с заданным именем
        string fileName = "staika_golubyov_step" + std::to_string(snap_number) + ".vtp";
        writeSpheresAsPoints(polyData, fileName);
    }
};

int main() {
    // кол-во голубей
    int N = 500;
    // Шаг по времени
    double tau = 0.01;
    // Время
    double T = 0;
    // Создаём сетку заданного размера
    CalcMesh mesh(N);

    // Пишем её начальное состояние в VTK
    mesh.snapshot(0);

    // Делаем шаги по времени, 
    // на каждом шаге считаем новое состояние и пишем его в VTK
    for (unsigned int step = 0; step < 5; step++) {
        for (int _step = 0; _step < 10; _step++){
            mesh.doTimeStep(tau, T);
	    T += tau;
        }
        
        mesh.snapshot(step);
    }

    return 0;
}
