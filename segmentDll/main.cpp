#define NULL 0
#include "../raytracing.h"
#include<list>
#include<vector>
#include<iostream>
#include <fstream>
#define DLL_API __declspec(dllexport)

using namespace std;
/*
int num_of_rays = 0;
bool beInit = false;
double* rays = NULL;
std::list<RAY_VOX*> segments;
std::vector<int> valid_ray_ids;
*/
long voxDimXY = 0;
long voxDimXYZ = 0;
int spDims[3] = { 0 };
int voxDims[3] = { 0 };
double *coords[3] = { NULL,NULL,NULL };
bool beSetVoxModel = false;
//ofstream ofile("d:/seg.txt");

/*
extern "C" DLL_API bool cmpfunction2(INDEX_LENGTH& a, INDEX_LENGTH & b) {
		return (a.id < b.id);
	}

extern "C"	DLL_API void segments_compuation_by_raytracing(double* rays, int numOfRay, double* split_surf_coords[3], int num_of_split_surfs[3], std::list<RAY_VOX*>& segments)
	{
		if (!beInit) return;
		valid_ray_ids.clear();
		segments.clear();
		//cout << "segments_compuation_by_raytracing0:" << endl;

		for (int i = 0; i < numOfRay; i++)
		{

			INDEX_LENGTH* newitc = NULL;
			if (i == 1795 || i==1080) {
				double minCoords[3] = { split_surf_coords[LONDIM][0],split_surf_coords[LATDIM][0],split_surf_coords[RDIM][0] };
				double maxCoords[3] = { split_surf_coords[LONDIM][num_of_split_surfs[LONDIM] - 1],split_surf_coords[LATDIM][num_of_split_surfs[LATDIM] - 1],split_surf_coords[RDIM][num_of_split_surfs[RDIM] - 1] };

				//cout << "Coords= " << rays[i * 6] << " " << rays[i * 6 + 1] << " " << rays[i * 6 + 2] << " " << rays[i * 6 + 3] << " " << rays[i * 6 + 4] << " " << rays[i * 6 + 5] << endl;
				//cout << minCoords[0] << " " << minCoords[1] << " " << minCoords[2] << " " << maxCoords[0] << " " << maxCoords[1] << " " << maxCoords[2] << endl;
			}

			int num = segments_of_a_ray_by_ray_tracing(rays + i * 6, rays + i * 6 + 3, split_surf_coords, num_of_split_surfs, newitc);

			if (num > 0) {
				RAY_VOX* new_ray_vox = new RAY_VOX();
				new_ray_vox->rayid = i;
				new_ray_vox->numOfVoxel = num;
				new_ray_vox->vox_itcs = newitc;
				segments.push_back(new_ray_vox);
				valid_ray_ids.push_back(i);
				//if (i == 1795 || i==1080)
				//	cout << ", " << new_ray_vox->rayid;

			}
			//cout << endl;
		}
	}

extern "C"	DLL_API void releaseRAY_VOXList() {
		for (auto e : segments) {
			if (e->vox_itcs) {
				delete[] e->vox_itcs;
				e->vox_itcs = NULL;
			}
		}
		segments.clear();
	}

extern "C"	DLL_API void releaseCoords() {
		if (coords[0]) {
			delete[] coords[0]; coords[0] = NULL;
		}
		if (coords[1]) {
			delete[] coords[1];
			coords[1] = NULL;
		}
		if (coords[2]) {
			delete[] coords[2];
			coords[2] = NULL;
		}

	}
*/
/*
extern "C"	DLL_API void destroy()
{
	releaseCoords();
	releaseRAY_VOXList();
	if (rays) delete[] rays;
	rays = NULL;
	valid_ray_ids.clear();
}
extern "C"	DLL_API void initialize(const double* lons, const double* lats, const double* rs, int nLon, int nLat, int nR, const double* xyz, int numOfRays)
	{
	//destroy();
		//cout << "numOfRays="<< numOfRays<<endl;
	destroy();
		num_of_rays = numOfRays;
		rays = new double[numOfRays * 6];
		for (int i = 0; i < numOfRays*6; i++) {
			rays[i] = xyz[i];
		}
		//cout << "nLon=" << nLon << "nLat=" << nLat << "nR=" << nR << endl;

		coords[0] = new double[nLon];
		coords[1] = new double[nLat];
		coords[2] = new double[nR];
		//cout << "coords[0]=" << endl;

		for (int i = 0; i < nLon; i++) {
			coords[0][i] = lons[i];
			//cout << "  " << coords[0][i] << endl;

		}
		//cout << "coords[1]=" << endl;

		for (int i = 0; i < nLat; i++) {
			coords[1][i] = lats[i];
			//cout << "  " << coords[1][i] << endl;

		}
		//cout << "coords[2]=" << endl;

		for (int i = 0; i < nR; i++) {
			coords[2][i] = rs[i];
			//cout << "  " << coords[2][i] << endl;

		}
		spDims[0] = nLon; spDims[1] = nLat; spDims[2] = nR;
		voxDims[0] = nLon - 1; voxDims[1] = nLat - 1; voxDims[2] = nR - 1;
		voxDimXY = (nLon - 1)*(nLat - 1);
		voxDimXYZ = voxDimXY * (nR - 1);
		//cout << "numOfRay=" << num_of_rays << endl;
		beInit = true;
	}
*/
extern "C"	DLL_API void release() {
	if (coords[0]) {
		delete[] coords[0]; coords[0] = NULL;
	}
	if (coords[1]) {
		delete[] coords[1];
		coords[1] = NULL;
	}
	if (coords[2]) {
		delete[] coords[2];
		coords[2] = NULL;
	}

}
extern "C"	DLL_API void setVoxelModel(const double* lons, const double* lats, const double* rs, int nLon, int nLat, int nR)
{
	//destroy();
		//cout << "numOfRays="<< numOfRays<<endl;
	release();


	coords[0] = new double[nLon];
	coords[1] = new double[nLat];
	coords[2] = new double[nR];
	//cout << "coords[0]=" << endl;

	for (int i = 0; i < nLon; i++) {
		coords[0][i] = lons[i];
		//cout << "  " << coords[0][i] << endl;

	}
	//cout << "coords[1]=" << endl;

	for (int i = 0; i < nLat; i++) {
		coords[1][i] = lats[i];
		//cout << "  " << coords[1][i] << endl;

	}
	//cout << "coords[2]=" << endl;

	for (int i = 0; i < nR; i++) {
		coords[2][i] = rs[i];
		//cout << "  " << coords[2][i] << endl;

	}
	spDims[0] = nLon; spDims[1] = nLat; spDims[2] = nR;
	voxDims[0] = nLon - 1; voxDims[1] = nLat - 1; voxDims[2] = nR - 1;
	voxDimXY = (nLon - 1)*(nLat - 1);
	voxDimXYZ = voxDimXY * (nR - 1);
	//cout << "nLon=" << nLon << "nLat=" << nLat << "nR=" << nR << endl;
	//cout << "spDims0=" << spDims[0] << "spDims1=" << spDims[1] << "spDims2=" << spDims[2] << endl;

	//cout << "numOfRay=" << num_of_rays << endl;
	beSetVoxModel = true;
}
/*
extern "C"	DLL_API int compute_segments()
	{
		if (!beInit) return -1;
		segments_compuation_by_raytracing(rays, num_of_rays,coords,spDims, segments);
		return segments.size();
	}

extern "C" DLL_API int get_valid_rays_num() {
	cout << "valid_ray_ids.size()=" << valid_ray_ids.size() << endl;
	if (valid_ray_ids.size()>0) return valid_ray_ids.size();
	double minCoords[3] = { coords[0][0],coords[1][0],coords[2][0] };
	double maxCoords[3] = { coords[0][spDims[0]-1],coords[1][spDims[1] - 1],coords[2][spDims[2] - 1] };
	valid_ray_ids.clear();
	cout << " get_valid_rays_num" <<",num_of_rays="<< num_of_rays<< endl;
	for (int i = 0; i < num_of_rays; i++) {
		//if (i == 1795 || i==1080) {
		//	cout << "\nCoords= " << rays[i*6]<< " " << rays[i * 6+1] << " " << rays[i * 6 + 2] << " " << rays[i * 6 + 3] << " " << rays[i * 6 + 4] << " " << rays[i * 6 + 5] << endl;
		//	cout << minCoords[0] << " " << minCoords[1] << " " << minCoords[2] << " " << maxCoords[0] << " " << maxCoords[1] << " " << maxCoords[2] << endl;
		//}
		if (is_a_valid_ray(rays + i * 6, rays + i * 6 + 3, minCoords, maxCoords, NULL))
		{
			cout << ", " << i;

			valid_ray_ids.push_back(i);
		}
	}
	cout << endl;

	return valid_ray_ids.size();
}
extern "C"	DLL_API int get_ray_ids(int* ids) {
	if (valid_ray_ids.size() < 1) {
		get_valid_rays_num();
	}
	int iCount = 0;

	for (std::vector<int>::iterator it = valid_ray_ids.begin(); it != valid_ray_ids.end(); it++)
	{
		ids[iCount++] = *it;
	}
	return valid_ray_ids.size();
}
extern "C"	DLL_API int get_segment_matrix(double* matrix)
	{
		if (segments.size() < 1){
			compute_segments();
		}
		cout << "segments.size()="<<segments.size()<<",voxDimXYZ=" << voxDimXYZ<< endl;
		//long iCount = 0;
		for (long j = 0; j < voxDimXYZ*segments.size(); j++) {
			matrix[j] = 0;
		}
		int iRay = 0;
		for (list<RAY_VOX*>::const_iterator it = segments.begin(); it != segments.end(); it++)
		{
			RAY_VOX* a_ray = *it;
			//sort(a_ray->vox_itcs, a_ray->vox_itcs + a_ray->numOfVoxel, cmpfunction2);
			//ofile << index;
			//long beginID = 0, endID = 0;
			//cout << "\nrayid:" << (*it)->rayid << ", voxNum" << a_ray->numOfVoxel << ",voxDimXYZ=" << voxDimXYZ << "endid=" << endID << endl;
			
			for (long j = 0; j < a_ray->numOfVoxel; j++) {

				if (a_ray->vox_itcs[j].id >= voxDimXYZ)
					cout << "voxID=" << a_ray->vox_itcs[j].id << ",voxDimXYZ=" << voxDimXYZ;
				matrix[iRay*voxDimXYZ + a_ray->vox_itcs[j].id] = a_ray->vox_itcs[j].intercept;

			}

			//for (long i = beginID; i < voxDimXYZ; i++)
			//{
			//	matrix[iCount++] = 0;
			//}
			iRay++;
		}
		return segments.size();
	}

*/
extern "C" DLL_API int get_segments_of_a_ray(double* xyzs,long long* voxid,double* segLength) {
	if (!beSetVoxModel) return 0;
	INDEX_LENGTH* newitc = NULL;
	long num = segments_of_a_ray_by_ray_tracing(xyzs, xyzs+ 3, coords, spDims, newitc);
	for (long i = 0; i < num; i++) {
		voxid[i] = newitc[i].index[2] * voxDimXY + newitc[i].index[1] * (spDims[0] - 1) + newitc[i].index[0];
		//ofile << voxid[i] << endl;
		segLength[i] = newitc[i].intercept;
	}
	return num;
}
extern "C" DLL_API int is_valid(double* xyzs) {
	if (!beSetVoxModel) return false;
	double minCoords[3] = { coords[0][0],coords[1][0],coords[2][0] };
	double maxCoords[3] = { coords[0][spDims[0] - 1],coords[1][spDims[1] - 1],coords[2][spDims[2] - 1] };
	if (is_a_valid_ray(xyzs, xyzs + 3, minCoords, maxCoords, NULL))
	{
		//cout << "True" << endl;
		return 1;
	}
	else {
		//cout << "False" << endl;

		return 0;
	}
}
extern "C" DLL_API void print_data() {
	cout << "voxDimXY=" << voxDimXY << endl;
	cout << "voxDimXYZ=" << voxDimXYZ << endl;
	cout << "spDims[3]=" << spDims[0]<<"," << spDims[1] << "," << spDims[2] << endl;
	cout << "voxDims[3]=" << voxDims[0] << "," << voxDims[1] << "," << voxDims[2] << endl;

	cout << "coords[0]=" << endl;
	if(coords[0])
	for (int i = 0; i < spDims[0]; i++) {
		cout << "  " << coords[0][i] << endl;

	}
	cout << "coords[1]=" << endl;
	if (coords[1])
	for (int i = 0; i < spDims[1]; i++) {
		cout << "  " << coords[1][i] << endl;
	}
	cout << "coords[2]=" << endl;
	if (coords[2])
	for (int i = 0; i < spDims[2]; i++) {
		//coords[2][i] = rs[i];
		cout << "  " << coords[2][i] << endl;
	}
	cout << "beSetVoxModel=" << beSetVoxModel;

}

