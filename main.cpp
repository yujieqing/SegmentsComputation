#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include "common.h"
#include "raytracing.h"
#include "tradition.h"
#include "pugixml.cpp"
#include <iomanip>
#include <stdlib.h>
#include <iostream>
#include<list>
#include <chrono>
#define NUMOFMODELS 10
using namespace std;

bool cmpfunction2(INDEX_LENGTH& a, INDEX_LENGTH & b) {
	return (a.id < b.id);
}
void read_rays_coords_from_file(const char* filename, double*& rays, int& numOfRay)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//function: read coordinates for rays from a files, and output the rays in an two-dimensional array, each row of the file records the two endpoints of the ray in longitude, latitude, and radial coordinates
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
{
	ifstream ifile(filename);
	vector<double> vals;
	string linebuf;
	while (getline(ifile, linebuf))
	{
		stringstream ss(linebuf);
		string tmp;
		while (getline(ss, tmp, ','))
		{
			vals.push_back(atof((tmp.c_str())));
		}
	
	}

	numOfRay = vals.size() / 6;

	if (rays) delete[] rays;
	rays = new double[numOfRay * 6];
	for (int i = 0; i < numOfRay; i++)
	{
		int j = i;// +begin;
		rays[i * 6] = vals[j * 6];
		rays[i * 6 + 1] = vals[j * 6 + 1];
		rays[i * 6 + 2] = vals[j * 6 + 2];
		rays[i * 6 + 3] = vals[j * 6 + 3];
		rays[i * 6 + 4] = vals[j * 6 + 4];
		rays[i * 6 + 5] = vals[j * 6 + 5];
	}
	ifile.close();

}
void read_voxel_model(const char* filename, double**& split_surf_coords, int num_of_split_surfs[3])
// read the information of  voxel model.
//      split_surf_coords: an two dimension array for the splitting surfaces' coordinates 
//      num_of_split_surfs: the number of splitting surfaces along each dimension
{

	double georadius;
	string voxeltype;
	pugi::xml_document doc;
	doc.load_file(filename);
	pugi::xml_node voxelmodel = doc.child("CITXML").child("VoxelModel");
	voxeltype = voxelmodel.attribute("type").as_string();
	bool usingdegree = false, usingHeight = false, usingkm = false;
	string angleUnit = voxelmodel.attribute("angleUnit").as_string();
	if (!strcmp(angleUnit.c_str(), "degree")) usingdegree = true;
	string distUnit = voxelmodel.attribute("distUnit").as_string();
	if (!strcmp(distUnit.c_str(), "km")) usingkm = true;
	string heightMode = voxelmodel.attribute("heightMode").as_string();;
	if (!strcmp(heightMode.c_str(), "height")) usingHeight = true;
	georadius = voxelmodel.attribute("radius").as_double();
	if (usingkm)		 georadius = georadius * 1000;
	split_surf_coords = new double*[3];
	if (!strcmp(voxeltype.c_str(), "UniformVoxel"))
	{
		double org[3], res[3];
		int dims[3];
		pugi::xml_node uniformmodel = voxelmodel.child("UniformVoxel");
		istringstream iss(uniformmodel.child("dimensions").first_child().value());
		iss >> dims[0] >> dims[1] >> dims[2];
		istringstream iss1(uniformmodel.child("origin").first_child().value());
		iss1 >> org[0] >> org[1] >> org[2];
		istringstream iss2(uniformmodel.child("resolution").first_child().value());
		iss2 >> res[0] >> res[1] >> res[2];

		if (usingdegree)
			for (int i = 0; i < 2; i++) {
				org[i] = org[i] / 180 * PI;
				res[i] = res[i] / 180 * PI;
			}
		if (usingkm)			 org[2] = org[2] * 1000;
		if (usingHeight) org[2] = org[2] + georadius;
		;
		num_of_split_surfs[LONDIM] = dims[LONDIM] + 1;
		num_of_split_surfs[LATDIM] = dims[LATDIM] + 1;
		num_of_split_surfs[RDIM] = dims[RDIM] + 1;
		split_surf_coords[LONDIM] = new double[num_of_split_surfs[LONDIM]];
		for (int i = 0; i < num_of_split_surfs[LONDIM]; i++)
			split_surf_coords[LONDIM][i] = org[0] + i * res[0];
		split_surf_coords[LATDIM] = new double[num_of_split_surfs[LATDIM]];
		for (int i = 0; i < num_of_split_surfs[LATDIM]; i++)
			split_surf_coords[LATDIM][i] = org[1] + i * res[1];
		split_surf_coords[RDIM] = new double[num_of_split_surfs[RDIM]];
		for (int i = 0; i < num_of_split_surfs[RDIM]; i++)
			split_surf_coords[RDIM][i] = org[2] + i * res[2];
	}
	else
	{
		pugi::xml_node rectmmodel = voxelmodel.child("RectilinearGrid");
		stringstream ss;
		double tmp;
		vector<double> lons, lats, rs;
		ss << rectmmodel.child("lonCoords").first_child().value();
		while (1) {
			ss >> tmp;
			if (ss.fail()) break;
			lons.push_back(tmp);
		}
		ss.clear();
		ss << rectmmodel.child("latCoords").first_child().value();
		while (1) {
			ss >> tmp;
			if (ss.fail()) break;
			lats.push_back(tmp);
		}
		ss.clear();
		ss << rectmmodel.child("rCoords").first_child().value();
		while (1) {
			ss >> tmp;
			if (ss.fail()) break;
			rs.push_back(tmp);
		}
		num_of_split_surfs[LONDIM] = lons.size();
		num_of_split_surfs[LATDIM] = lats.size();
		num_of_split_surfs[RDIM] = rs.size();
		split_surf_coords[LONDIM] = new double[num_of_split_surfs[LONDIM]];
		split_surf_coords[LATDIM] = new double[num_of_split_surfs[LATDIM]];
		split_surf_coords[RDIM] = new double[num_of_split_surfs[RDIM]];
		for (int i = 0; i < num_of_split_surfs[LONDIM]; i++) {
			split_surf_coords[LONDIM][i] = lons[i];
			if (usingdegree)	 split_surf_coords[LONDIM][i] = split_surf_coords[LONDIM][i] / 180 * PI;

		}

		for (int i = 0; i < num_of_split_surfs[LATDIM]; i++) {
			double lat = lats[i];
			split_surf_coords[LATDIM][i] = lat;
			if (usingdegree)	split_surf_coords[LATDIM][i] = split_surf_coords[LATDIM][i] / 180 * PI;
		}

		for (int i = 0; i < num_of_split_surfs[RDIM]; i++) {
			split_surf_coords[RDIM][i] = rs[i];
			if (usingkm)		 split_surf_coords[RDIM][i] *= 1000;
			if (usingHeight)	 split_surf_coords[RDIM][i] += georadius;
		}
	}

}
void segments_compuation_by_raytracing(double* rays, int numOfRay, double* split_surf_coords[3], int num_of_split_surfs[3], std::list<RAY_VOX*>& segments)
//compute segements for a number of rays using ray traing method
{
	for (int i = 0; i < numOfRay; i++)
	{
		INDEX_LENGTH* newitc = NULL;
		int num = segments_of_a_ray_by_ray_tracing(rays + i * 6, rays + i * 6 + 3, split_surf_coords, num_of_split_surfs, newitc);
		if (num > 0) {
			RAY_VOX* new_ray_vox = new RAY_VOX();
			new_ray_vox->rayid = i;
			new_ray_vox->numOfVoxel = num;
			new_ray_vox->vox_itcs = newitc;
			segments.push_back(new_ray_vox);
		}
	}
}

void segments_computation_by_tradition(double* rays, int numOfRay, double* split_surf_coords[3], int num_of_split_surfs[3], std::list<RAY_VOX*>& segments)
//compute segements for a number of rays using ray traing method

{
	for (int i = 0; i < numOfRay; i++)
	{
		vector<INDEX_LENGTH*> segments_of_a_ray;
		int num = segments_of_a_ray_by_tradition(rays + i * 6, rays + i * 6 + 3, split_surf_coords, num_of_split_surfs, segments_of_a_ray);
		if (num > 0) {
			RAY_VOX* new_ray_vox = new RAY_VOX();
			new_ray_vox->rayid = i;
			new_ray_vox->numOfVoxel = num;
			INDEX_LENGTH* segments_of_a_ray_Ptr = new INDEX_LENGTH[segments_of_a_ray.size()];
			new_ray_vox->vox_itcs = segments_of_a_ray_Ptr;
			// the following codes can be removed when comparing with ray tracing method. They just copy data to a new structure
			for (int j = 0; j < segments_of_a_ray.size(); j++)
			{
				segments_of_a_ray_Ptr[j].intercept = segments_of_a_ray[j]->intercept;
				segments_of_a_ray_Ptr[j].t_vec_first_pt = segments_of_a_ray[j]->t_vec_first_pt;
				segments_of_a_ray_Ptr[j].index[0] = segments_of_a_ray[j]->index[0];
				segments_of_a_ray_Ptr[j].index[1] = segments_of_a_ray[j]->index[1];
				segments_of_a_ray_Ptr[j].index[2] = segments_of_a_ray[j]->index[2];
				segments_of_a_ray_Ptr[j].id = segments_of_a_ray[j]->id;
			}			
			segments.push_back(new_ray_vox);
		}

	}
}

void save_segments_by_matrix(const char* filename, std::list<RAY_VOX*>& segments,long voxDims)
// save the resutant segments in form of matix 
{
	ofstream ofile(filename);
	for (list<RAY_VOX*>::const_iterator it = segments.begin(); it != segments.end(); it++)
	{
		RAY_VOX* a_ray = *it;
		sort(a_ray->vox_itcs, a_ray->vox_itcs + a_ray->numOfVoxel, cmpfunction2);
		//ofile << index;
		long beginID = 0, endID = 0;
		for (long j = 0; j<a_ray->numOfVoxel; j++) {

			endID = a_ray->vox_itcs[j].id;
			for (int i = beginID; i < endID; i++)
			{
				if (i == voxDims - 1)
					ofile << 0;
				else
					ofile << 0 << ",";
			}
			if (endID == voxDims - 1)
			ofile << a_ray->vox_itcs[j].intercept;
			else
				ofile << a_ray->vox_itcs[j].intercept << ",";

			beginID = endID+1;
		}
		for (int i = beginID; i < voxDims; i++)
		{
			if (i == voxDims - 1)
				ofile << 0;
			else
				ofile << 0 << ",";
		}
		ofile << endl;
	}
	ofile.close();

}

void save_segments_by_cent_size_val(const char* filename, std::list<RAY_VOX*>& segments, double** split_surf_coords, int num_of_split_surfs[3])
// save the segments in center, size and length format.Each row of the file is a segment in a voxel. There are seven columns in a row
// The first six columns indicate the center and size of a voxel by spherical coordinates(in radians)
// the last column indicates the length of the segment associated with the voxel
{
	int numOfRay = segments.size();
	int dimVox10 = (num_of_split_surfs[LATDIM] - 1)*(num_of_split_surfs[LONDIM] - 1);
	ofstream ofile(filename);
	for (list<RAY_VOX*>::const_iterator it = segments.begin(); it != segments.end(); it++)
	{

		RAY_VOX* a_ray = *it;
		sort(a_ray->vox_itcs, a_ray->vox_itcs + a_ray->numOfVoxel, cmpfunction2);

		for (int j = 0; j < a_ray->numOfVoxel; j++)
		{
			long id = a_ray->vox_itcs[j].id;
			long rIndex = id / dimVox10;
			long restIndex = id % dimVox10;
			long latIndex = restIndex / (num_of_split_surfs[LONDIM] - 1);
			long lonIndex = restIndex % (num_of_split_surfs[LONDIM] - 1);
			ofile << (split_surf_coords[LONDIM][lonIndex + 1] + split_surf_coords[LONDIM][lonIndex]) / 2 << ", "
				<< (split_surf_coords[LATDIM][latIndex + 1] + split_surf_coords[LATDIM][latIndex]) / 2 << ", "
				<< (split_surf_coords[RDIM][rIndex + 1] + split_surf_coords[RDIM][rIndex]) / 2 << ", "
				<< (split_surf_coords[LONDIM][lonIndex + 1] - split_surf_coords[LONDIM][lonIndex]) << ", "
				<< (split_surf_coords[LATDIM][latIndex + 1] - split_surf_coords[LATDIM][latIndex]) << ", "
				<< (split_surf_coords[RDIM][rIndex + 1] - split_surf_coords[RDIM][rIndex]) << ", "
				<< a_ray->vox_itcs[j].intercept	<< endl;
		}
	}
	ofile.close();
}
void save_segments_by_sparse_matrix(const char* filename, std::list<RAY_VOX*>& segments)
// Each row of the file is a segment in a voxel. There are three colomns in a row, i.e., ray ID, voxel ID, segment's length
{
	int numOfRay = segments.size();
	ofstream ofile(filename);
	for (list<RAY_VOX*>::const_iterator it = segments.begin(); it != segments.end(); it++)
	{
		RAY_VOX* a_ray = *it;
		sort(a_ray->vox_itcs, a_ray->vox_itcs + a_ray->numOfVoxel, cmpfunction2);
		for (int j = 0; j < a_ray->numOfVoxel; j++)
		{
			long id = a_ray->vox_itcs[j].id;
			ofile << a_ray->rayid<<","<< id << ", "<< a_ray->vox_itcs[j].intercept << endl;
		}
	}
	ofile.close();
}
void save_segments_by_indice(const char* filename, std::list<RAY_VOX*>& segments, double** split_surf_coords, int num_of_split_surfs[3])
// ray ID, voxel ID, longitudinal index, latitudinal index, radial index,  segment's length
{
	int numOfRay = segments.size();
	int dimVox10 = (num_of_split_surfs[LATDIM] - 1)*(num_of_split_surfs[LONDIM] - 1);
	ofstream ofile(filename);
	for (list<RAY_VOX*>::const_iterator it = segments.begin(); it != segments.end(); it++)
	{

		RAY_VOX* a_ray = *it;
		sort(a_ray->vox_itcs, a_ray->vox_itcs + a_ray->numOfVoxel, cmpfunction2);
		for (int j = 0; j < a_ray->numOfVoxel; j++)
		{
			long id = a_ray->vox_itcs[j].id;
			
			long rIndex = id / dimVox10;
			long restIndex = id % dimVox10;
			long latIndex = restIndex / (num_of_split_surfs[LONDIM] - 1);
			long lonIndex = restIndex % (num_of_split_surfs[LONDIM] - 1);
			ofile << a_ray->rayid << ","
				<< id << ", "
				<< lonIndex << ", "
				<< latIndex << ", "
				<< rIndex << ", "
				<< a_ray->vox_itcs[j].intercept << endl;
		}
	}
	ofile.close();
}

int main()
{
	string ray_file = "input/rayCoordinates.txt";
	string path = "input/voxelModels/";
	string voxelFiles[NUMOFMODELS] = {
		"1.xml","2.xml","3.xml","4.xml","5.xml","6.xml","7.xml","8.xml","9.xml","10.xml"
	};
	const char runtimeFile[]= "output/time_cost.txt";
	ofstream ofile(runtimeFile);
	ofile << "id, raytracing,tradition" << endl;
	double* rays = NULL;
	int numOfRay=0;
	read_rays_coords_from_file(ray_file.c_str(), rays, numOfRay);
	for (int i = 0; i < NUMOFMODELS; i++)
	{
		ofile << i << ", ";
		std::list<RAY_VOX*> segments;
		double ** split_surf_coords = NULL;
		int num_of_split_surfs[3];
		read_voxel_model((path+voxelFiles[i]).c_str(), split_surf_coords, num_of_split_surfs);
		std::cout << "start ray tracing job " << i<<endl;
		chrono::high_resolution_clock::time_point tt1 = chrono::high_resolution_clock::now();
		segments_compuation_by_raytracing(rays, numOfRay, split_surf_coords, num_of_split_surfs, segments);
		chrono::high_resolution_clock::time_point tt2 = chrono::high_resolution_clock::now();

		auto time_span1 = chrono::duration_cast<chrono::microseconds>(tt2 - tt1);
		double ms = time_span1.count() / 1000.0;
		std::cout << "\tray tracing cost " <<ms<< "ms" <<endl;
		ofile <<  ms<<",";
		char filename[1024];
		sprintf_s(filename, 1024, "output/raytracing_centSize%d.txt", i);
		long voxdim = (num_of_split_surfs[LONDIM] - 1)* (num_of_split_surfs[LATDIM] - 1)* (num_of_split_surfs[RDIM] - 1);
		save_segments_by_cent_size_val(filename, segments, split_surf_coords, num_of_split_surfs);
		sprintf_s(filename, 1024, "output/raytracing_sparse%d.txt", i);
		save_segments_by_sparse_matrix(filename, segments);

		std::cout << "start tradition job " <<i<< endl;
		std::list<RAY_VOX*> segments1;

		tt1 = chrono::high_resolution_clock::now();
		segments_computation_by_tradition(rays, numOfRay, split_surf_coords, num_of_split_surfs, segments1);
		tt2 = chrono::high_resolution_clock::now();
		auto time_span2 = chrono::duration_cast<chrono::microseconds>(tt2 - tt1);
		ms = time_span2.count() / 1000.0;
		ofile << ms << endl;
		std::cout << "\ttradition cost " << ms << "ms" << endl;
		sprintf_s(filename, 1024, "output/tradition_centSize%d.txt", i);
		save_segments_by_cent_size_val(filename, segments1, split_surf_coords, num_of_split_surfs);
		sprintf_s(filename, 1024, "output/tradition_sparse%d.txt", i);
		save_segments_by_sparse_matrix(filename, segments1);

	}
	ofile.close();
}