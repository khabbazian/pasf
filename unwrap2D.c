/*2D phase unwrapping, modified for inclusion in scipy by Gregor Thalhammer

This program was written by Munther Gdeisat and Miguel Arevallilo Herraez to program the two-dimensional unwrapper
entitled "Fast two-dimensional phase-unwrapping algorithm based on sorting by
reliability following a noncontinuous path"
by  Miguel Arevallilo Herraez, David R. Burton, Michael J. Lalor, and Munther A. Gdeisat
published in the Journal Applied Optics, Vol. 41, No. 35, pp. 7437, 2002.
This program was written by Munther Gdeisat, Liverpool John Moores University, United Kingdom.
Date 26th August 2007
The wrapped phase map is assumed to be of floating point data type. The resultant unwrapped phase map is also of floating point type.
The mask is of byte data type.
When the mask is 255 this means that the pixel is valid
When the mask is 0 this means that the pixel is invalid (noisy or corrupted pixel)
This program takes into consideration the image wrap around problem encountered in MRI imaging.*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mex.h"


static double PI = 3.141592654f;
static double TWOPI = 6.283185307f;

#define NOMASK 0
#define MASK 1

typedef struct
{
	double mod;
	int x_connectivity;
	int y_connectivity;
	int no_of_edges;
} params_t;

struct PIXELM
{
	int increment;		        
	int number_of_pixels_in_group;  
	double value;			
	double reliability;
	bool input_mask;	
	bool extended_mask;	
	int group;			
	int new_group;
	struct PIXELM *head;		
	struct PIXELM *last;		
	struct PIXELM *next;		
};

typedef struct PIXELM PIXELM;

struct EDGE
{
	double reliab;		
	PIXELM *pointer_1;	
	PIXELM *pointer_2;	
	int increment;		
};

typedef struct EDGE EDGE;

#define swap(x,y) {EDGE t; t=x; x=y; y=t;}
#define order(x,y) if (x.reliab > y.reliab) swap(x,y)
#define o2(x,y) order(x,y)
#define o3(x,y,z) o2(x,y); o2(x,z); o2(y,z)

typedef enum {yes, no} yes_no;

yes_no find_pivot(EDGE *left, EDGE *right, double *pivot_ptr)
{
	EDGE a, b, c, *p;

	a = *left;
	b = *(left + (right - left) /2 );
	c = *right;
	o3(a,b,c);

	if (a.reliab < b.reliab)
	{
		*pivot_ptr = b.reliab;
		return yes;
	}

	if (b.reliab < c.reliab)
	{
		*pivot_ptr = c.reliab;
		return yes;
	}

	for (p = left + 1; p <= right; ++p)
	{
		if (p->reliab != left->reliab)
		{
			*pivot_ptr = (p->reliab < left->reliab) ? left->reliab : p->reliab;
			return yes;
		}
		return no;
	}

	return no;
}

EDGE *partition(EDGE *left, EDGE *right, double pivot)
{
	while (left <= right)
	{
		while (left->reliab < pivot)
			++left;
		while (right->reliab >= pivot)
			--right;
		if (left < right)
		{
			swap (*left, *right);
			++left;
			--right;
		}
	}
	return left;
}

void quicker_sort(EDGE *left, EDGE *right)
{
	EDGE *p;
	double pivot;

	if (find_pivot(left, right, &pivot) == yes)
	{
		p = partition(left, right, pivot);
		quicker_sort(left, p - 1);
		quicker_sort(p, right);
	}
}

void  initialisePIXELs(double *wrapped_image, bool *input_mask, bool *extended_mask, PIXELM *pixel, int image_width, int image_height)
{
	PIXELM *pixel_pointer = pixel;
	double *wrapped_image_pointer = wrapped_image;
	bool *input_mask_pointer = input_mask;
	bool *extended_mask_pointer = extended_mask;
	int i, j;

	for (i=0; i < image_height; i++)
	{
		for (j=0; j < image_width; j++)
		{
			pixel_pointer->increment = 0;
			pixel_pointer->number_of_pixels_in_group = 1;
			pixel_pointer->value = *wrapped_image_pointer;
			pixel_pointer->reliability = 9999999.f + rand();
			pixel_pointer->input_mask = *input_mask_pointer;
			pixel_pointer->extended_mask = *extended_mask_pointer;
			pixel_pointer->head = pixel_pointer;
			pixel_pointer->last = pixel_pointer;
			pixel_pointer->next = NULL;
			pixel_pointer->new_group = 0;
			pixel_pointer->group = -1;
			pixel_pointer++;
			wrapped_image_pointer++;
			input_mask_pointer++;
			extended_mask_pointer++;
		}
	}
}

double wrap(double pixel_value)
{
	double wrapped_pixel_value;
	if (pixel_value > PI)	wrapped_pixel_value = pixel_value - TWOPI;
	else if (pixel_value < -PI) wrapped_pixel_value = pixel_value + TWOPI;
	else wrapped_pixel_value = pixel_value;
	return wrapped_pixel_value;
}


int find_wrap(double pixelL_value, double pixelR_value)
{
	double difference;
	int wrap_value;
	difference = pixelL_value - pixelR_value;

	if (difference > PI)	wrap_value = -1;
	else if (difference < -PI)	wrap_value = 1;
	else wrap_value = 0;

	return wrap_value;
}

void extend_mask(bool *input_mask, bool *extended_mask,
		int image_width, int image_height,
		params_t *params)
{
	int i,j;
	int image_width_plus_one = image_width + 1;
	int image_width_minus_one = image_width - 1;
	bool *IMP = input_mask    + image_width + 1;	
	bool *EMP = extended_mask + image_width + 1;	

	for (i=1; i < image_height - 1; ++i)
	{
		for (j=1; j < image_width - 1; ++j)
		{
			if ( (*IMP) == NOMASK && (*(IMP + 1) == NOMASK) && (*(IMP - 1) == NOMASK) &&
					(*(IMP + image_width) == NOMASK) && (*(IMP - image_width) == NOMASK) &&
					(*(IMP - image_width_minus_one) == NOMASK) && (*(IMP - image_width_plus_one) == NOMASK) &&
					(*(IMP + image_width_minus_one) == NOMASK) && (*(IMP + image_width_plus_one) == NOMASK) )
			{
				*EMP = NOMASK;
			}
			++EMP;
			++IMP;
		}
		EMP += 2;
		IMP += 2;
	}

	if (params->x_connectivity == 1)
	{
		IMP = input_mask    + 2 * image_width - 1;
		EMP = extended_mask + 2 * image_width -1;
		for (i=1; i < image_height - 1; ++ i)
		{
			if ( (*IMP) == NOMASK && (*(IMP - 1) == NOMASK) &&  (*(IMP + 1) == NOMASK) &&
					(*(IMP + image_width) == NOMASK) && (*(IMP - image_width) == NOMASK) &&
					(*(IMP - image_width - 1) == NOMASK) && (*(IMP - image_width + 1) == NOMASK) &&
					(*(IMP + image_width - 1) == NOMASK) && (*(IMP - 2 * image_width + 1) == NOMASK) )
			{
				*EMP = NOMASK;
			}
			EMP += image_width;
			IMP += image_width;
		}

		IMP = input_mask    + image_width;
		EMP = extended_mask + image_width;
		for (i=1; i < image_height - 1; ++i)
		{
			if ( (*IMP) == NOMASK && (*(IMP - 1) == NOMASK) && (*(IMP + 1) == NOMASK) &&
					(*(IMP + image_width) == NOMASK) && (*(IMP - image_width) == NOMASK) &&
					(*(IMP - image_width + 1) == NOMASK) && (*(IMP + image_width + 1) == NOMASK) &&
					(*(IMP + image_width - 1) == NOMASK) && (*(IMP + 2 * image_width - 1) == NOMASK) )
			{
				*EMP = NOMASK;
			}
			EMP += image_width;
			IMP += image_width;
		}
	}

	if (params->y_connectivity == 1)
	{
		IMP = input_mask    + 1;
		EMP = extended_mask + 1;
		for (i=1; i < image_width - 1; ++i)
		{
			if ( (*IMP) == NOMASK && (*(IMP - 1) == NOMASK) && (*(IMP + 1) == NOMASK) &&
					(*(IMP + image_width) == NOMASK) && (*(IMP + image_width * (image_height - 1)) == NOMASK) &&
					(*(IMP + image_width + 1) == NOMASK) && (*(IMP + image_width - 1) == NOMASK) &&
					(*(IMP + image_width * (image_height - 1) - 1) == NOMASK) && (*(IMP + image_width * (image_height - 1) + 1) == NOMASK) )
			{
				*EMP = NOMASK;
			}
			EMP++;
			IMP++;
		}

		IMP = input_mask    + image_width * (image_height - 1) + 1;
		EMP = extended_mask + image_width * (image_height - 1) + 1;
		for (i=1; i < image_width - 1; ++i)
		{
			if ( (*IMP) == NOMASK && (*(IMP - 1) == NOMASK) && (*(IMP + 1) == NOMASK) &&
					(*(IMP - image_width) == NOMASK) && (*(IMP - image_width - 1) == NOMASK) && (*(IMP - image_width + 1) == NOMASK) &&
					(*(IMP - image_width * (image_height - 1)    ) == NOMASK) &&
					(*(IMP - image_width * (image_height - 1) - 1) == NOMASK) &&
					(*(IMP - image_width * (image_height - 1) + 1) == NOMASK) )
			{
				*EMP = NOMASK;
			}
			EMP++;
			IMP++;
		}
	}
}

void calculate_reliability(double *wrappedImage, PIXELM *pixel,
		int image_width, int image_height,
		params_t *params)
{
	int image_width_plus_one = image_width + 1;
	int image_width_minus_one = image_width - 1;
	PIXELM *pixel_pointer = pixel + image_width_plus_one;
	double *WIP = wrappedImage + image_width_plus_one; 
	double H, V, D1, D2;
	int i, j;

	for (i = 1; i < image_height -1; ++i)
	{
		for (j = 1; j < image_width - 1; ++j)
		{
			if (pixel_pointer->extended_mask == NOMASK)
			{
				H = wrap(*(WIP - 1) - *WIP) - wrap(*WIP - *(WIP + 1));
				V = wrap(*(WIP - image_width) - *WIP) - wrap(*WIP - *(WIP + image_width));
				D1 = wrap(*(WIP - image_width_plus_one) - *WIP) - wrap(*WIP - *(WIP + image_width_plus_one));
				D2 = wrap(*(WIP - image_width_minus_one) - *WIP) - wrap(*WIP - *(WIP + image_width_minus_one));
				pixel_pointer->reliability = H*H + V*V + D1*D1 + D2*D2;
			}
			pixel_pointer++;
			WIP++;
		}
		pixel_pointer += 2;
		WIP += 2;
	}

	if (params->x_connectivity == 1)
	{
		PIXELM *pixel_pointer = pixel + image_width;
		double *WIP = wrappedImage + image_width;

		for (i = 1; i < image_height - 1; ++i)
		{
			if (pixel_pointer->extended_mask == NOMASK)
			{
				H = wrap(*(WIP + image_width - 1) - *WIP) - wrap(*WIP - *(WIP + 1));
				V = wrap(*(WIP - image_width) - *WIP) - wrap(*WIP - *(WIP + image_width));
				D1 = wrap(*(WIP - 1) - *WIP) - wrap(*WIP - *(WIP + image_width_plus_one));
				D2 = wrap(*(WIP - image_width_minus_one) - *WIP) - wrap(*WIP - *(WIP + 2* image_width - 1));
				pixel_pointer->reliability = H*H + V*V + D1*D1 + D2*D2;
			}
			pixel_pointer += image_width;
			WIP += image_width;
		}

		pixel_pointer = pixel + 2 * image_width - 1;
		WIP = wrappedImage + 2 * image_width - 1;

		for (i = 1; i < image_height - 1; ++i)
		{
			if (pixel_pointer->extended_mask == NOMASK)
			{
				H = wrap(*(WIP - 1) - *WIP) - wrap(*WIP - *(WIP - image_width_minus_one));
				V = wrap(*(WIP - image_width) - *WIP) - wrap(*WIP - *(WIP + image_width));
				D1 = wrap(*(WIP - image_width_plus_one) - *WIP) - wrap(*WIP - *(WIP + 1));
				D2 = wrap(*(WIP - 2 * image_width - 1) - *WIP) - wrap(*WIP - *(WIP + image_width_minus_one));
				pixel_pointer->reliability = H*H + V*V + D1*D1 + D2*D2;
			}
			pixel_pointer += image_width;
			WIP += image_width;
		}
	}

	if (params->y_connectivity == 1)
	{
		PIXELM *pixel_pointer = pixel + 1;
		double *WIP = wrappedImage + 1;

		for (i = 1; i < image_width - 1; ++i)
		{
			if (pixel_pointer->extended_mask == NOMASK)
			{
				H =  wrap(*(WIP - 1) - *WIP) - wrap(*WIP - *(WIP + 1));
				V =  wrap(*(WIP + image_width*(image_height - 1)) - *WIP) - wrap(*WIP - *(WIP + image_width));
				D1 = wrap(*(WIP + image_width*(image_height - 1) - 1) - *WIP) - wrap(*WIP - *(WIP + image_width_plus_one));
				D2 = wrap(*(WIP + image_width*(image_height - 1) + 1) - *WIP) - wrap(*WIP - *(WIP + image_width_minus_one));
				pixel_pointer->reliability = H*H + V*V + D1*D1 + D2*D2;
			}
			pixel_pointer++;
			WIP++;
		}

		pixel_pointer = pixel + (image_height - 1) * image_width + 1;
		WIP = wrappedImage + (image_height - 1) * image_width + 1;

		for (i = 1; i < image_width - 1; ++i)
		{
			if (pixel_pointer->extended_mask == NOMASK)
			{
				H =  wrap(*(WIP - 1) - *WIP) - wrap(*WIP - *(WIP + 1));
				V =  wrap(*(WIP - image_width) - *WIP) - wrap(*WIP - *(WIP -(image_height - 1) * (image_width)));
				D1 = wrap(*(WIP - image_width_plus_one) - *WIP) - wrap(*WIP - *(WIP - (image_height - 1) * (image_width) + 1));
				D2 = wrap(*(WIP - image_width_minus_one) - *WIP) - wrap(*WIP - *(WIP - (image_height - 1) * (image_width) - 1));
				pixel_pointer->reliability = H*H + V*V + D1*D1 + D2*D2;
			}
			pixel_pointer++;
			WIP++;
		}
	}
}

void  horizontalEDGEs(PIXELM *pixel, EDGE *edge,
		int image_width, int image_height,
		params_t *params)
{
	int i, j;
	EDGE *edge_pointer = edge;
	PIXELM *pixel_pointer = pixel;
	int no_of_edges = params->no_of_edges;

	for (i = 0; i < image_height; i++)
	{
		for (j = 0; j < image_width - 1; j++)
		{
			if (pixel_pointer->input_mask == NOMASK && (pixel_pointer + 1)->input_mask == NOMASK)
			{
				edge_pointer->pointer_1 = pixel_pointer;
				edge_pointer->pointer_2 = (pixel_pointer+1);
				edge_pointer->reliab = pixel_pointer->reliability + (pixel_pointer + 1)->reliability;
				edge_pointer->increment = find_wrap(pixel_pointer->value, (pixel_pointer + 1)->value);
				edge_pointer++;
				no_of_edges++;
			}
			pixel_pointer++;
		}
		pixel_pointer++;
	}
	if (params->x_connectivity == 1)
	{
		pixel_pointer = pixel + image_width - 1;
		for (i = 0; i < image_height; i++)
		{
			if (pixel_pointer->input_mask == NOMASK && (pixel_pointer - image_width + 1)->input_mask == NOMASK)
			{
				edge_pointer->pointer_1 = pixel_pointer;
				edge_pointer->pointer_2 = (pixel_pointer - image_width + 1);
				edge_pointer->reliab = pixel_pointer->reliability + (pixel_pointer - image_width + 1)->reliability;
				edge_pointer->increment = find_wrap(pixel_pointer->value, (pixel_pointer  - image_width + 1)->value);
				edge_pointer++;
				no_of_edges++;
			}
			pixel_pointer+=image_width;
		}
	}
	params->no_of_edges = no_of_edges;
}

void  verticalEDGEs(PIXELM *pixel, EDGE *edge,
		int image_width, int image_height,
		params_t *params)
{
	int i, j;
	int no_of_edges = params->no_of_edges;
	PIXELM *pixel_pointer = pixel;
	EDGE *edge_pointer = edge + no_of_edges;

	for (i=0; i < image_height - 1; i++)
	{
		for (j=0; j < image_width; j++)
		{
			if (pixel_pointer->input_mask == NOMASK && (pixel_pointer + image_width)->input_mask == NOMASK)
			{
				edge_pointer->pointer_1 = pixel_pointer;
				edge_pointer->pointer_2 = (pixel_pointer + image_width);
				edge_pointer->reliab = pixel_pointer->reliability + (pixel_pointer + image_width)->reliability;
				edge_pointer->increment = find_wrap(pixel_pointer->value, (pixel_pointer + image_width)->value);
				edge_pointer++;
				no_of_edges++;
			}
			pixel_pointer++;
		} 
	} 

	if (params->y_connectivity == 1)
	{
		pixel_pointer = pixel + image_width *(image_height - 1);
		for (i = 0; i < image_width; i++)
		{
			if (pixel_pointer->input_mask == NOMASK && (pixel_pointer - image_width *(image_height - 1))->input_mask == NOMASK)
			{
				edge_pointer->pointer_1 = pixel_pointer;
				edge_pointer->pointer_2 = (pixel_pointer - image_width *(image_height - 1));
				edge_pointer->reliab = pixel_pointer->reliability + (pixel_pointer - image_width *(image_height - 1))->reliability;
				edge_pointer->increment = find_wrap(pixel_pointer->value, (pixel_pointer - image_width *(image_height - 1))->value);
				edge_pointer++;
				no_of_edges++;
			}
			pixel_pointer++;
		}
	}
	params->no_of_edges = no_of_edges;
}

void  gatherPIXELs(EDGE *edge, params_t *params)
{
	int k;
	PIXELM *PIXEL1;
	PIXELM *PIXEL2;
	PIXELM *group1;
	PIXELM *group2;
	EDGE *pointer_edge = edge;
	int incremento;

	for (k = 0; k < params->no_of_edges; k++)
	{
		PIXEL1 = pointer_edge->pointer_1;
		PIXEL2 = pointer_edge->pointer_2;

		if (PIXEL2->head != PIXEL1->head)
		{
			if ((PIXEL2->next == NULL) && (PIXEL2->head == PIXEL2))
			{
				PIXEL1->head->last->next = PIXEL2;
				PIXEL1->head->last = PIXEL2;
				(PIXEL1->head->number_of_pixels_in_group)++;
				PIXEL2->head=PIXEL1->head;
				PIXEL2->increment = PIXEL1->increment-pointer_edge->increment;
			}

			else if ((PIXEL1->next == NULL) && (PIXEL1->head == PIXEL1))
			{
				PIXEL2->head->last->next = PIXEL1;
				PIXEL2->head->last = PIXEL1;
				(PIXEL2->head->number_of_pixels_in_group)++;
				PIXEL1->head = PIXEL2->head;
				PIXEL1->increment = PIXEL2->increment+pointer_edge->increment;
			}

			else
			{
				group1 = PIXEL1->head;
				group2 = PIXEL2->head;
				if (group1->number_of_pixels_in_group > group2->number_of_pixels_in_group)
				{
					group1->last->next = group2;
					group1->last = group2->last;
					group1->number_of_pixels_in_group = group1->number_of_pixels_in_group + group2->number_of_pixels_in_group;
					incremento = PIXEL1->increment-pointer_edge->increment - PIXEL2->increment;
					while (group2 != NULL)
					{
						group2->head = group1;
						group2->increment += incremento;
						group2 = group2->next;
					}
				}

				else
				{
					group2->last->next = group1;
					group2->last = group1->last;
					group2->number_of_pixels_in_group = group2->number_of_pixels_in_group + group1->number_of_pixels_in_group;
					incremento = PIXEL2->increment + pointer_edge->increment - PIXEL1->increment;
					while (group1 != NULL)
					{
						group1->head = group2;
						group1->increment += incremento;
						group1 = group1->next;
					}

				}
			}
		}
		pointer_edge++;
	}
}

void  unwrapImage(PIXELM *pixel, int image_width, int image_height)
{
	int i;
	int image_size = image_width * image_height;
	PIXELM *pixel_pointer=pixel;

	for (i = 0; i < image_size; i++)
	{
		pixel_pointer->value += TWOPI * (double)(pixel_pointer->increment);
		pixel_pointer++;
	}
}

void  maskImage(PIXELM *pixel, bool *input_mask, int image_width, int image_height)
{
	PIXELM *pointer_pixel = pixel;
	bool *IMP = input_mask;	
	double min=99999999.f;
	int i;
	int image_size = image_width * image_height;

	for (i = 0; i < image_size; i++)
	{
		if ((pointer_pixel->value < min) && (*IMP == NOMASK))
			min = pointer_pixel->value;

		pointer_pixel++;
		IMP++;
	}

	pointer_pixel = pixel;
	IMP = input_mask;

	for (i = 0; i < image_size; i++)
	{
		if ((*IMP) == MASK)
		{
			pointer_pixel->value = min;
		}
		pointer_pixel++;
		IMP++;
	}
}

void  returnImage(PIXELM *pixel, double *unwrapped_image, int image_width, int image_height)
{
	int i;
	int image_size = image_width * image_height;
	double *unwrapped_image_pointer = unwrapped_image;
	PIXELM *pixel_pointer = pixel;

	for (i=0; i < image_size; i++)
	{
		*unwrapped_image_pointer = pixel_pointer->value;
		pixel_pointer++;
		unwrapped_image_pointer++;
	}
}

	void
unwrap2D(double* wrapped_image, double* UnwrappedImage, bool* input_mask,
		int image_width, int image_height,
		int wrap_around_x, int wrap_around_y)
{
	params_t params = {TWOPI, wrap_around_x, wrap_around_y, 0};
	bool *extended_mask;
	PIXELM *pixel;
	EDGE *edge;
	int image_size = image_height * image_width;
	int No_of_Edges_initially = 2 * image_width * image_height;

	extended_mask = (bool *) calloc(image_size, sizeof(bool));
	pixel = (PIXELM *) calloc(image_size, sizeof(PIXELM));
	edge = (EDGE *) calloc(No_of_Edges_initially, sizeof(EDGE));

	extend_mask(input_mask, extended_mask, image_width, image_height, &params);
	initialisePIXELs(wrapped_image, input_mask, extended_mask, pixel, image_width, image_height);
	calculate_reliability(wrapped_image, pixel, image_width, image_height, &params);
	horizontalEDGEs(pixel, edge, image_width, image_height, &params);
	verticalEDGEs(pixel, edge, image_width, image_height, &params);

	quicker_sort(edge, edge + params.no_of_edges - 1);

	gatherPIXELs(edge, &params);

	unwrapImage(pixel, image_width, image_height);
	maskImage(pixel, input_mask, image_width, image_height);

	returnImage(pixel, UnwrappedImage, image_width, image_height);

	free(edge);
	free(pixel);
	free(extended_mask);
}


void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
	double sum; 	
	/* Check for proper number of arguments. */
	if(nrhs != 2) {
		mexErrMsgTxt("Two inputs required.");
	} else if(nlhs>1) {
		mexErrMsgTxt("Too many output arguments.");
	}

	/*Find number of rows and columns in input data*/
	int rows = mxGetM(prhs[0]);		
	int cols = mxGetN(prhs[0]);

	/*Get pointer to input data*/
	double* wrapped = mxGetPr(prhs[0]);
	bool* mask  = (bool*) mxGetPr(prhs[1]);

	if( mxGetM(prhs[1]) != rows || mxGetN(prhs[1]) != cols )
		mexErrMsgTxt("Phase and mask array have to be of the same size.");

	/*int i=0;
	for(i=0; i<rows*cols; ++i)
		mask[i]=0;*/

	/*create space for output*/
	plhs[0] = mxCreateDoubleMatrix(rows, cols, mxREAL);

	/*Get pointer to output array*/
	double* unwrapped = mxGetPr(plhs[0]);

	/*unwrap2D(wrapped, unwrapped, mask, rows, cols, false, false); */
	unwrap2D(wrapped, unwrapped, mask, cols, rows, false, false);
}
