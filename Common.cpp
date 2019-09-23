#pragma once
#include"common.h"
#include <math.h>
int t_vec_to_sphere(const double pt1[3], const double pt2[3], const double rCoord, double t_vector[2])
{
	double a = pt2[XDIM] - pt1[XDIM];
	double b = pt2[YDIM] - pt1[YDIM];
	double c = pt2[ZDIM] - pt1[ZDIM];
	double A = a * a + b * b + c * c;
	double B = 2 * (a*pt1[XDIM] + b * pt1[YDIM] + c * pt1[ZDIM]);
	double C = pt1[XDIM] * pt1[XDIM] + pt1[YDIM] * pt1[YDIM] + pt1[ZDIM] * pt1[ZDIM] - rCoord * rCoord;
	double bac = B * B - 4 * A*C;
	if (bac < 0) return 0;
	double deta = sqrt(B*B - 4 * A*C);
	double t[2];
	t[0] = (-B - deta) / 2 / A;
	t[1] = (-B + deta) / 2 / A;
	int index = 0;
	if (t[0] >= -THRESHOLD_T_OF_TWO_POINTS && t[0] <= 1 + THRESHOLD_T_OF_TWO_POINTS)
		t_vector[index++] = t[0];
	if (t[1] >= -THRESHOLD_T_OF_TWO_POINTS && t[1] <= 1 + THRESHOLD_T_OF_TWO_POINTS)
		t_vector[index++] = t[1];
	return index;
}
int t_vec_to_latitude_surface(const double pt1[3], const double pt2[3], const double latCoord, double t_vector[2])
//latCoord: in radian
{
	double a = pt2[XDIM] - pt1[XDIM];
	double b = pt2[YDIM] - pt1[YDIM];
	double c = pt2[ZDIM] - pt1[ZDIM];
	double k = tan(latCoord);
	k = k * k;
	double A = k * (a*a + b * b) - c * c;
	double B = 2 * (k*(a*pt1[XDIM] + b * pt1[YDIM]) - c * pt1[ZDIM]);
	double C = k * pt1[XDIM] * pt1[XDIM] + k * pt1[YDIM] * pt1[YDIM] - pt1[ZDIM] * pt1[ZDIM];
	double bac = B * B - 4 * A*C;
	if (bac < 0) return 0;
	double deta = sqrt(bac);
	double t[2];
	t[0] = (-B - deta) / 2 / A;
	t[1] = (-B + deta) / 2 / A;
	int index = 0;
	if (t[0] >= -THRESHOLD_T_OF_TWO_POINTS && t[0] <= 1 + THRESHOLD_T_OF_TWO_POINTS) t_vector[index++] = t[0];
	if (t[1] >= -THRESHOLD_T_OF_TWO_POINTS && t[1] <= 1 + THRESHOLD_T_OF_TWO_POINTS) t_vector[index++] = t[1];
	return index;
}
int t_vec_to_longitude_surface(const double pt1[3], const double pt2[3], const double lonCoord, double& t_vector)
//lonCoord:in radian
{
	double a = pt2[XDIM] - pt1[XDIM];
	double b = pt2[YDIM] - pt1[YDIM];
	double c = pt2[ZDIM] - pt1[ZDIM];
	double cosa = cos(lonCoord);
	double sina = sin(lonCoord);
	double t = (pt1[YDIM] * cosa - pt1[XDIM] * sina) / (a*sina - b * cosa);
	if (t >= -THRESHOLD_T_OF_TWO_POINTS && t <= 1 + THRESHOLD_T_OF_TWO_POINTS) { t_vector = t; return 1; }
	return 0;
}
void XYZ_from_t_vec(const double pt1[3], const double pt2[3], const double t, double result[3])
{
	result[XDIM] = pt1[XDIM] + t * (pt2[XDIM] - pt1[XDIM]);
	result[YDIM] = pt1[YDIM] + t * (pt2[YDIM] - pt1[YDIM]);
	result[ZDIM] = pt1[ZDIM] + t * (pt2[ZDIM] - pt1[ZDIM]);
}
void  LLR_to_XYZ(const double oldCoord[3], double newCoord[3])
{
	newCoord[ZDIM] = oldCoord[RDIM] * sin(oldCoord[LATDIM]);
	newCoord[YDIM] = oldCoord[RDIM] * cos(oldCoord[LATDIM])*sin(oldCoord[LONDIM]);
	newCoord[XDIM] = oldCoord[RDIM] * cos(oldCoord[LATDIM])*cos(oldCoord[LONDIM]);
}

void XYZ_to_LLR(const double oldCoord[3], double newCoord[3])
{
	newCoord[RDIM] = sqrt(oldCoord[XDIM] * oldCoord[XDIM] + oldCoord[YDIM] * oldCoord[YDIM] + oldCoord[ZDIM] * oldCoord[ZDIM]);
	newCoord[LATDIM] = atan2(oldCoord[ZDIM], sqrt(oldCoord[XDIM] * oldCoord[XDIM] + oldCoord[YDIM] * oldCoord[YDIM]));
	newCoord[LONDIM] = atan2(oldCoord[YDIM], oldCoord[XDIM]);
}
bool is_between_lon_range(double lower_boundary, double upper_boundary, double lon)
//is a longitude value between two given longitude values
{
	double lon_10 = upper_boundary - lower_boundary;
	if (lon_10 < 0) lon_10 += 2 * PI;
	if (lon_10 >= PI + PI) lon_10 = lon_10 - PI - PI;
	double dlon = lon - lower_boundary;
	if (dlon < 0) dlon += 2 * PI;
	if (dlon >= PI + PI) dlon = dlon - PI - PI;

	return (dlon <= lon_10);
}

bool is_a_valid_ray(double pt1[3], double pt2[3], double minCoords[3], double maxCoords[3], double* t_firstLast)
{
	/////////////////////////////////////////////////////////////////////////////
	//1�����������붥������Ľ��㼰��Ӧ��t vector
	//2�������޳���Ч����
	/////////////////////////////////////////////////////////////////////////////
	double tt_low[2], tt_high[2], lt[2], low_pt_t, high_pt_t, low_pt_xyz[3], high_pt_xyz[3], high_pt_llr[3], low_pt_llr[3];
	int num_t_low_sphere = t_vec_to_sphere(pt1, pt2, minCoords[RDIM], tt_low);
	int num_t_high_sphere = t_vec_to_sphere(pt1, pt2, maxCoords[RDIM], tt_high);


	//�ᴩ��������
	if (num_t_low_sphere == 1 && num_t_high_sphere == 1) { low_pt_t = tt_low[0];	high_pt_t = tt_high[0]; }
	//���ڲ����򶥲����棬���ɶ��������ڲ�
	else if (num_t_low_sphere == 0 && num_t_high_sphere == 1) {
		double r1 = sqrt(pt1[0] * pt1[0] + pt1[1] * pt1[1] + pt1[2] * pt1[2]);
		if (r1 > minCoords[RDIM] && r1 < maxCoords[RDIM])
			//p1λ�ڴ���������
		{
			low_pt_t = 0;
			high_pt_t = tt_high[0];
		}
		else//p2λ�ڴ���������
		{
			low_pt_t = 1;
			high_pt_t = tt_low[0];
		}
	}
	//���ڲ�����ײ����棬���ɵײ������ڲ� 
	else if (num_t_low_sphere == 1 && num_t_high_sphere == 0) {
		double r1 = sqrt(pt1[0] * pt1[0] + pt1[1] * pt1[1] + pt1[2] * pt1[2]);
		if (r1 > minCoords[RDIM] && r1 < maxCoords[RDIM])
			//��p1λ�ڴ���������
		{
			low_pt_t = tt_low[0];
			high_pt_t = 0;
		}
		else//p2λ�ڴ���������
		{
			low_pt_t = tt_low[0];
			high_pt_t = 1;
		}
	}
	//ˮƽ������������
	else if (num_t_low_sphere == 0 && num_t_high_sphere == 2) { low_pt_t = tt_high[0];	high_pt_t = tt_high[1]; }
	//�����Ϊ��Ч����
	else	return false;

	XYZ_from_t_vec(pt1, pt2, low_pt_t, low_pt_xyz);
	XYZ_from_t_vec(pt1, pt2, high_pt_t, high_pt_xyz);
	XYZ_to_LLR(high_pt_xyz, high_pt_llr);
	XYZ_to_LLR(low_pt_xyz, low_pt_llr);

	//��������ɸѡ�󣬻������½�һ��ɸѡ
	//1)���������������о�����������ͶӰ��Χ֮�⣬ҲΪ��Ч����
	if (!(is_between_lon_range(minCoords[LONDIM], maxCoords[LONDIM], high_pt_llr[LONDIM])
		&& is_between_lon_range(minCoords[LONDIM], maxCoords[LONDIM], low_pt_llr[LONDIM])
		&& high_pt_llr[LATDIM] >= minCoords[LATDIM] && high_pt_llr[LATDIM] <= maxCoords[LATDIM]
		&& low_pt_llr[LATDIM] >= minCoords[LATDIM] && low_pt_llr[LATDIM] <= maxCoords[LATDIM]
		)) return false;

	//2)ֻҪ�������ϱ�γ��������κ�2�����㼴��Ч���ߣ���ҪΪ���ų�������ˮƽ����ʼ��γ����ͬ������ĳһγ���棬��Ϊγ���������棬������ֱ��
	if (t_vec_to_latitude_surface(pt1, pt2, minCoords[LATDIM], lt) == 2) return false;
	if (t_vec_to_latitude_surface(pt1, pt2, maxCoords[LATDIM], lt) == 2) return false;

	if (t_firstLast)
	{
		t_firstLast[0] = low_pt_t;
		t_firstLast[1] = high_pt_t;
	}
	return true;
	/////////////////////////////////////////////////////////////////////////////
}