#include "tradition.h"
#include <assert.h>
#include <algorithm>
#include <vector>

bool cmpfunction1(double& a, double & b) {
	return (a < b);
}

double segment_length_in_a_voxel(double minCoords[3], double maxCoords[3], double startPt[3], double endPt[3], double vectNormal)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//function: calucate the length of a ray inside a voxel
//Input:
//		minCoords: the minimum value of each dimension for a voxel
//      maxCoords: the maximum value of each dimension for a voxel
//		startPt: the start point of a ray
//      endPt: the end point of a ray
//		vectNormal: the length of the ray's vector
//Output:
//      orientVector: the length of the ray inside the voxel

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
{
	double t_vals[10];
	int nCount = 0;

	double tmp_t[2], tmp_XYZ[3], tmp_LLR[3];

	int tmp_num;
	// intersetion with the first longitudinal surface 
	if (t_vec_to_longitude_surface(startPt, endPt, minCoords[LONDIM], tmp_t[0]) > 0)
	{
		XYZ_from_t_vec(startPt, endPt, tmp_t[0], tmp_XYZ);
		XYZ_to_LLR(tmp_XYZ, tmp_LLR);
		if (tmp_LLR[LATDIM] >= minCoords[LATDIM] && tmp_LLR[LATDIM] <= maxCoords[LATDIM]
			&& tmp_LLR[RDIM] >= minCoords[RDIM] && tmp_LLR[RDIM] <= maxCoords[RDIM])
		{
			t_vals[nCount] = tmp_t[0];
			nCount++;
		}
	}
	// intersetion with the second longitudinal surface 
	if (t_vec_to_longitude_surface(startPt, endPt, maxCoords[LONDIM], tmp_t[0]) > 0)
	{
		XYZ_from_t_vec(startPt, endPt, tmp_t[0], tmp_XYZ);
		XYZ_to_LLR(tmp_XYZ, tmp_LLR);
		if (tmp_LLR[LATDIM] >= minCoords[LATDIM] && tmp_LLR[LATDIM] <= maxCoords[LATDIM]
			&& tmp_LLR[RDIM] >= minCoords[RDIM] && tmp_LLR[RDIM] <= maxCoords[RDIM])
		{
			t_vals[nCount] = tmp_t[0];
			nCount++;
		}
	}

	// intersetion with the first latitudinal surface 
	tmp_num = t_vec_to_latitude_surface(startPt, endPt, minCoords[LATDIM], tmp_t);
	for (int i = 0; i < tmp_num; i++)
	{
		XYZ_from_t_vec(startPt, endPt, tmp_t[i], tmp_XYZ);
		XYZ_to_LLR(tmp_XYZ, tmp_LLR);
		if (tmp_LLR[RDIM] >= minCoords[RDIM] && tmp_LLR[RDIM] <= maxCoords[RDIM]
			&& is_between_lon_range(minCoords[LONDIM], maxCoords[LONDIM], tmp_LLR[LONDIM]))
		{
			t_vals[nCount] = tmp_t[i];
			nCount++;
		}
	}

	// intersetion with the second latitudinal surface 
	tmp_num = t_vec_to_latitude_surface(startPt, endPt, maxCoords[LATDIM], tmp_t);
	for (int i = 0; i < tmp_num; i++)
	{
		XYZ_from_t_vec(startPt, endPt, tmp_t[i], tmp_XYZ);
		XYZ_to_LLR(tmp_XYZ, tmp_LLR);
		if (tmp_LLR[RDIM] >= minCoords[RDIM] && tmp_LLR[RDIM] <= maxCoords[RDIM]
			&& is_between_lon_range(minCoords[LONDIM], maxCoords[LONDIM], tmp_LLR[LONDIM]))
		{
			t_vals[nCount] = tmp_t[i];
			nCount++;
		}
	}

	// intersetion with the first spherical surface 
	tmp_num = t_vec_to_sphere(startPt, endPt, minCoords[RDIM], tmp_t);
	for (int i = 0; i < tmp_num; i++)
	{
		XYZ_from_t_vec(startPt, endPt, tmp_t[i], tmp_XYZ);
		XYZ_to_LLR(tmp_XYZ, tmp_LLR);
		if (tmp_LLR[LATDIM] >= minCoords[LATDIM] && tmp_LLR[LATDIM] <= maxCoords[LATDIM]
			&& is_between_lon_range(minCoords[LONDIM], maxCoords[LONDIM], tmp_LLR[LONDIM]))
		{
			t_vals[nCount] = tmp_t[i];
			nCount++;
		}
	}

	// intersetion with the second spherical surface 
	tmp_num = t_vec_to_sphere(startPt, endPt, maxCoords[RDIM], tmp_t);
	for (int i = 0; i < tmp_num; i++)
	{
		XYZ_from_t_vec(startPt, endPt, tmp_t[i], tmp_XYZ);
		XYZ_to_LLR(tmp_XYZ, tmp_LLR);
		if (tmp_LLR[LATDIM] >= minCoords[LATDIM] && tmp_LLR[LATDIM] <= maxCoords[LATDIM]
			&& is_between_lon_range(minCoords[LONDIM], maxCoords[LONDIM], tmp_LLR[LONDIM]))
		{
			t_vals[nCount] = tmp_t[i];
			nCount++;
		}
	}
	if (nCount < 2) return 0;
	if (nCount > 6) return 0;//there are six possible intersections maximally. Under this case, the two intersections are two vertexs of the intersecting voxel themself
	if (nCount == 2) 
		return fabs(t_vals[1] - t_vals[0])*vectNormal;
	else //merge those intersections that are very close to each other
	{
		std::sort(t_vals, t_vals + nCount, cmpfunction1);
		double t_finals[6];
		int iCurr_final = 0;
		t_finals[iCurr_final] = t_vals[0];
		for (int i = 1; i < nCount; i++)
		{
			if (fabs(t_finals[iCurr_final] - t_vals[i]) > THRESHOLD_T_OF_TWO_POINTS)
			{
				iCurr_final++;
				t_finals[iCurr_final] = t_vals[i];
			}
			i++;
		}
		assert(iCurr_final > 1);// If this is happened, we should increase the THRESHOLD_T_OF_TWO_POINTS
		return fabs(t_finals[1] - t_finals[0])*vectNormal;
	}
}
int segments_of_a_ray_by_tradition(double startPt[3], double endPt[3], double* split_surf_coords[3], int num_of_split_surfs[3], std::vector<INDEX_LENGTH*>& segments)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//function: calucate the segments' length using tradition method
//Input:
//		startPt: start point of a ray
//      endPt: end point of a ray
//      split_surf_coords: an two dimension array for the splitting surfaces' coordinates 
//      num_of_split_surfs: the number of splitting surfaces along each dimension
//Output:
//      segments: the segments' length using tradition method

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
{
	double minCoords[3] = { split_surf_coords[LONDIM][0],split_surf_coords[LATDIM][0],split_surf_coords[RDIM][0] };
	double maxCoords[3] = { split_surf_coords[LONDIM][num_of_split_surfs[LONDIM] - 1],split_surf_coords[LATDIM][num_of_split_surfs[LATDIM] - 1],split_surf_coords[RDIM][num_of_split_surfs[RDIM] - 1] };
	if (!is_a_valid_ray(startPt, endPt, minCoords, maxCoords, NULL)) return 0;


	double normal = sqrt((startPt[XDIM] - endPt[XDIM])*(startPt[XDIM] - endPt[XDIM]) + (startPt[YDIM] - endPt[YDIM])
		*(startPt[YDIM] - endPt[YDIM]) + (startPt[ZDIM] - endPt[ZDIM])*(startPt[ZDIM] - endPt[ZDIM]));
	double len = 0;
	long nCount = 0;
	long dim10 = (num_of_split_surfs[LATDIM] - 1) * (num_of_split_surfs[LONDIM] - 1);

	for (int k = 0; k < num_of_split_surfs[RDIM] - 1; k++)
	{
		minCoords[RDIM] = split_surf_coords[RDIM][k];
		maxCoords[RDIM] = split_surf_coords[RDIM][k + 1];
		for (int j = 0; j < num_of_split_surfs[LATDIM] - 1; j++)
		{
			minCoords[LATDIM] = split_surf_coords[LATDIM][j];
			maxCoords[LATDIM] = split_surf_coords[LATDIM][j + 1];
			for (int i = 0; i < num_of_split_surfs[LONDIM] - 1; i++)
			{
				minCoords[LONDIM] = split_surf_coords[LONDIM][i];
				maxCoords[LONDIM] = split_surf_coords[LONDIM][i + 1];
				if (i == 7 and j == 4 and k == 10)
					int a = 0;
				if ((len = segment_length_in_a_voxel(minCoords, maxCoords, startPt, endPt, normal)) > 0)
				{
					INDEX_LENGTH* pnew = new INDEX_LENGTH;
					pnew->intercept = len;
					pnew->index[LONDIM] = i;
					pnew->index[LATDIM] = j;
					pnew->index[RDIM] = k;
					pnew->id = pnew->index[RDIM] * dim10 + pnew->index[LATDIM] * (num_of_split_surfs[LONDIM] - 1) + pnew->index[LONDIM];
					segments.push_back(pnew);
					nCount++;
				}
			}
		}
	}
	return nCount;
}
