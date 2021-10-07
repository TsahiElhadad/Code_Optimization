//	206214165	Tzahi Elhadad
#include <stdbool.h> 

#define MIN(a,b) ((a)<(b)?(a):(b)) // macro methods more efficent.

typedef struct {
	unsigned char red;
	unsigned char green;
	unsigned char blue;
} pixel;

/* blur - function for blur the pictrue.
* Apply the kernel over each pixel.
* Ignore pixels where the kernel exceeds bounds. These are pixels with row index smaller than kernelSize/2 and/or
* column index smaller than kernelSize/2
*/
void blur(int kernelScale, pixel* src, Image* image) {
	// The keyword "register" hints the compiler that a given variable can be put in a register,
	// make the locallaty increased.

	int dim, kernelSize; // for not pass 7 arguments.
	dim = m; // for not pass 7 arguments.
	kernelSize = 3; // for not pass 7 arguments.
	int end, start;
	start = 1; // start line of pixcel check
	end = dim - start; // end line of pixcel check
	register int red = 0, green = 0, blue = 0; // sum of coloros to calculate
	register pixel* dst = (pixel*)image->data; // pointer to image data and convert to pixcel
	register pixel the_pixel, p1, p2, p3; // pointers to columns in row to check
	register int cal_indexes_range = 0; // position in 2d array that implemented in 1d array.
	int i, j, ii;
	register int my_index = 0; // position of index we check
	register int mult_out_from_for1 = 0, mult_out_from_for2 = 0;
	for (i = start; i < end; i++) {

		ii = i - 1; // calculate out of loop more efficent
		mult_out_from_for1 = (ii) * (dim); // calculate out of loop more efficent
		mult_out_from_for2 = i * dim; // calculate out of loop more efficent

		for (j = start; j < end; j++) {
			red = 0, green = 0, blue = 0;

			my_index = mult_out_from_for2 + j; // calculate the position of the current pixcel

			// Calculate all the array values around me for a 3 by 3 size matrix in one loop

			cal_indexes_range = mult_out_from_for1 +(j - 1); // calculate direction of first column in 3*3 mat
			p1 = src[cal_indexes_range]; // first column
			p2 = src[cal_indexes_range + 1]; // second column
			p3 = src[cal_indexes_range + 2]; // third column
			// calculate the first line colors 
			red += ((int)p1.red) + ((int)p2.red) + ((int)p3.red); // Instead of call sum_pixels_by_weight
			green += ((int)p1.green) + ((int)p2.green) + ((int)p3.green); // Instead of call sum_pixels_by_weight
			blue += ((int)p1.blue) + ((int)p2.blue) + ((int)p3.blue); // Instead of call sum_pixels_by_weight

			cal_indexes_range += dim; // move next line
			p1 = src[cal_indexes_range]; // first column
			p2 = src[cal_indexes_range + 1]; // second column
			p3 = src[cal_indexes_range + 2]; // third column

			// calculate the second line colors 
			red += ((int)p1.red) + ((int)p2.red) + ((int)p3.red); // Instead of call sum_pixels_by_weight
			green += ((int)p1.green) + ((int)p2.green) + ((int)p3.green); // Instead of call sum_pixels_by_weight
			blue += ((int)p1.blue) + ((int)p2.blue) + ((int)p3.blue); // Instead of call sum_pixels_by_weight

			cal_indexes_range += dim; // move next line
			p1 = src[cal_indexes_range]; // first column
			p2 = src[cal_indexes_range + 1]; // second column
			p3 = src[cal_indexes_range + 2]; // third column

			// calculate the last line colors 
			red += ((int)p1.red) + ((int)p2.red) + ((int)p3.red); // Instead of call sum_pixels_by_weight
			green += ((int)p1.green) + ((int)p2.green) + ((int)p3.green); // Instead of call sum_pixels_by_weight
			blue += ((int)p1.blue) + ((int)p2.blue) + ((int)p3.blue); // Instead of call sum_pixels_by_weight

			the_pixel.red = red / kernelScale;
			the_pixel.blue = blue / kernelScale;
			the_pixel.green = green / kernelScale;

			dst[my_index] = the_pixel;
		}
	}
}

// Function that sharpen the image.
void sharpen(int kernelScale, pixel* src, Image* image) {
	// The keyword register hints to compiler that a given variable can be put in a register
	// make the locallaty increased.

	int dim =m, kernelSize = 3; // m global variable and allways 3 for kernel size
	int end, start; // variables for run over the image matrix
	start = 1; // we can't calculate pixcels in the corners
	end = dim - start; // end of the seace index in the image matrix
	register int red = 0, green = 0, blue = 0; // for sum the colors in the 3*3 matrix
	register pixel* dst = (pixel*)image->data; // pointer to image data and convert to pixcel

	register pixel the_pixel; // the pixcel we calculate
	register pixel pixel_1, pixel_2, pixel_3; // pointers to columns in row to check.
	register int cal_indexes_range = 0; // start column of matrix 3*3 --> (i-1,j-1)
	int i, j, ii;
	register int my_index = 0, mult_line_dim = 0, mult_ii_dim = 0;
	for (i = start; i < end; i++) {

		ii = i - 1; // sum out of loop more efficent.
		mult_line_dim = i * dim; // multiple out of loop more efficent.
		mult_ii_dim = (ii) * (dim); // multiple out of loop more efficent.

		for (j = start; j < end; j++) {
			// Calculate all the array values around me for a 3 by 3 size matrix in one loop

			red = 0, green = 0, blue = 0; // reset sum
			my_index = mult_line_dim + j;  // calculate the pixcel location we are

			cal_indexes_range = mult_ii_dim +(j - 1); // start column of matrix 3*3 --> (i-1,j-1)
			pixel_1 = src[cal_indexes_range]; // pointer to first col
			pixel_2 = src[cal_indexes_range + 1]; // pointer to second col
			pixel_3 = src[cal_indexes_range + 2]; // pointer to last col
			// sum of colors
			red += ((int)pixel_1.red * -1) + ((int)pixel_2.red * -1) + ((int)pixel_3.red * -1); // Instead of call sum_pixels_by_weight
			green += ((int)pixel_1.green * -1) + ((int)pixel_2.green * -1) + ((int)pixel_3.green * -1); // Instead of call sum_pixels_by_weight
			blue += ((int)pixel_1.blue * -1) + ((int)pixel_2.blue * -1) + ((int)pixel_3.blue * -1); // Instead of call sum_pixels_by_weight

			cal_indexes_range += dim; // move next raw
			pixel_1 = src[cal_indexes_range];
			pixel_2 = src[cal_indexes_range + 1];
			pixel_3 = src[cal_indexes_range + 2];
			// sum of colors
			red += ((int)pixel_1.red * -1) + ((int)pixel_2.red * 9) + ((int)pixel_3.red * -1); // Instead of call sum_pixels_by_weight
			green += ((int)pixel_1.green * -1) + ((int)pixel_2.green * 9) + ((int)pixel_3.green * -1); // Instead of call sum_pixels_by_weight
			blue += ((int)pixel_1.blue * -1) + ((int)pixel_2.blue * 9) + ((int)pixel_3.blue * -1); // Instead of call sum_pixels_by_weight

			cal_indexes_range += dim; // move next raw and last
			pixel_1 = src[cal_indexes_range];
			pixel_2 = src[cal_indexes_range + 1];
			pixel_3 = src[cal_indexes_range + 2];
			// sum of colors
			red += ((int)pixel_1.red * -1) + ((int)pixel_2.red * -1) + ((int)pixel_3.red * -1); // Instead of call sum_pixels_by_weight
			green += ((int)pixel_1.green * -1) + ((int)pixel_2.green * -1) + ((int)pixel_3.green * -1); // Instead of call sum_pixels_by_weight
			blue += ((int)pixel_1.blue * -1) + ((int)pixel_2.blue * -1) + ((int)pixel_3.blue * -1); // Instead of call sum_pixels_by_weight
			
			the_pixel.red = (unsigned char)(red < 0 ? 0 : MIN(255, red));// only one call to min
			the_pixel.green = (unsigned char)(green < 0 ? 0 : MIN(255, green));// only one call to min
			the_pixel.blue = (unsigned char)(blue < 0 ? 0 : MIN(255, blue));// only one call to min
			
			dst[my_index] = the_pixel;
		}
	}
}

// Function for high filtered picture, run when bool filter is on(1).
void blurFilter(int kernelScale, pixel* src, Image* image) {
	// The keyword register hints to compiler that a given variable can be put in a register
	// make the locallaty increased.

	int dim = m, kernelSize = 3; // m global variable and allways 3 for kernel size
	int end, start; // variables for run over the image matrix
	start = 1; // we can't calculate pixcels in the corners
	end = dim - start; // we can't calculate pixcels in the corners
	register int red = 0, green = 0, blue = 0; // for sum the colors in the 3*3 matrix
	register pixel* dst = (pixel*)image->data; // pointer to image data and convert to pixcel
	register pixel the_pixel; // the current pixcel we calculate.
	register pixel pix_1, pix_2, pix_3; // pointers to columns pixcels in row to check.
	register int cal_indexes_range = 0; // start column of matrix 3*3 --> (i-1,j-1)
	int i, j, ii;
	int my_index = 0; // location of the pixel we calculate
	register int sum1 = 0, sum2 = 0, sum3 = 0; // sum of each pixcel[i,j] in the row.
	int min_p = 0, max_p = 0; // for calculate the maximum and minimum indexes(row, col) of 3*3 matrix check.
	register int tempMax = -1; // arbitrary value that is lower than minimum possible intensity, which is 0
	register int tempMin = 766; // arbitrary value that is higher than maximum possible intensity, which is 255*3=765
	register int mult_i_dim = 0, mult_ii_dim = 0;
	for (i = start; i < end; i++) {

		ii = i - 1; // multiple only here out of for
		mult_i_dim = i * dim; // multiple only here out of for
		mult_ii_dim = (ii) * (dim); // multiple only here out of for

		for (j = start; j < end; j++) {
			// Calculate all the array values around me for a 3 by 3 size matrix in one loop
			// calculate each row for increase locallaty and more calculates in one itreation.

			tempMax = -1, tempMin = 766; // reset the min and max each loop
			red = 0, green = 0, blue = 0; // reset the min and max each loop
			my_index = mult_i_dim + j; // calculate the pixcel location we are

			cal_indexes_range = mult_ii_dim +(j - 1); // calculate direction of first column in 3*3 mat
			pix_1 = src[cal_indexes_range]; // first column - miss
			pix_2 = src[cal_indexes_range + 1]; // second column - hit
			pix_3 = src[cal_indexes_range + 2]; // third column - hit
			// calculate the sum of each pixcel in the raw
			sum1 = ((int)pix_1.red) + ((int)pix_1.blue) + ((int)pix_1.green);
			sum2 = ((int)pix_2.red) + ((int)pix_2.blue) + ((int)pix_2.green);
			sum3 = ((int)pix_3.red) + ((int)pix_3.blue) + ((int)pix_3.green);
			// Get Minimum sum values of pixcel
			if (sum1 <= tempMin) {
				tempMin = sum1;
				min_p = cal_indexes_range;
			}
			if (sum2 <= tempMin) {
				tempMin = sum2;
				min_p = cal_indexes_range + 1;
			}
			if (sum3 <= tempMin) {
				tempMin = sum3;
				min_p = cal_indexes_range + 2;
			} // Get Maximum sum values of pixcel
			if (sum1 > tempMax) {
				tempMax = sum1;
				max_p = cal_indexes_range;
			}
			if (sum2 > tempMax) {
				tempMax = sum2;
				max_p = cal_indexes_range + 1;
			}
			if (sum3 > tempMax) {
				tempMax = sum3;
				max_p = cal_indexes_range + 2;
			}
			// sum of colors
			red += ((int)pix_1.red) + ((int)pix_2.red) + ((int)pix_3.red); // Instead of call sum_pixels_by_weight
			green += ((int)pix_1.green) + ((int)pix_2.green) + ((int)pix_3.green); // Instead of call sum_pixels_by_weight
			blue += ((int)pix_1.blue) + ((int)pix_2.blue) + ((int)pix_3.blue); // Instead of call sum_pixels_by_weight

			cal_indexes_range += dim; // move next line
			pix_1 = src[cal_indexes_range];
			pix_2 = src[cal_indexes_range + 1];
			pix_3 = src[cal_indexes_range + 2];
			// calculate the sum of each pixcel in the raw
			sum1 = ((int)pix_1.red) + ((int)pix_1.blue) + ((int)pix_1.green);
			sum2 = ((int)pix_2.red) + ((int)pix_2.blue) + ((int)pix_2.green);
			sum3 = ((int)pix_3.red) + ((int)pix_3.blue) + ((int)pix_3.green);
			// Get Minimum sum values of pixcel
			if (sum1 <= tempMin) {
				tempMin = sum1;
				min_p = cal_indexes_range;
			}
			if (sum2 <= tempMin) {
				tempMin = sum2;
				min_p = cal_indexes_range + 1;
			}
			if (sum3 <= tempMin) {
				tempMin = sum3;
				min_p = cal_indexes_range + 2;
			} // Get Maximum sum values of pixcel
			if (sum1 > tempMax) {
				tempMax = sum1;
				max_p = cal_indexes_range;
			}
			if (sum2 > tempMax) {
				tempMax = sum2;
				max_p = cal_indexes_range + 1;
			}
			if (sum3 > tempMax) {
				tempMax = sum3;
				max_p = cal_indexes_range + 2;
			}
			// sum of colors
			red += ((int)pix_1.red) + ((int)pix_2.red) + ((int)pix_3.red); // Instead of call sum_pixels_by_weight
			green += ((int)pix_1.green) + ((int)pix_2.green) + ((int)pix_3.green); // Instead of call sum_pixels_by_weight
			blue += ((int)pix_1.blue) + ((int)pix_2.blue) + ((int)pix_3.blue); // Instead of call sum_pixels_by_weight

			cal_indexes_range += dim; // move next line
			pix_1 = src[cal_indexes_range];
			pix_2 = src[cal_indexes_range + 1];
			pix_3 = src[cal_indexes_range + 2];
			// calculate the sum of each pixcel in the raw
			sum1 = ((int)pix_1.red) + ((int)pix_1.blue) + ((int)pix_1.green);
			sum2 = ((int)pix_2.red) + ((int)pix_2.blue) + ((int)pix_2.green);
			sum3 = ((int)pix_3.red) + ((int)pix_3.blue) + ((int)pix_3.green);
			// Get Minimum sum values of pixcel
			if (sum1 <= tempMin) {
				tempMin = sum1;
				min_p = cal_indexes_range;
			}
			if (sum2 <= tempMin) {
				tempMin = sum2;
				min_p = cal_indexes_range + 1;
			}
			if (sum3 <= tempMin) {
				tempMin = sum3;
				min_p = cal_indexes_range + 2;
			} // Get Maximum sum values of pixcel
			if (sum1 > tempMax) {
				tempMax = sum1;
				max_p = cal_indexes_range;
			}
			if (sum2 > tempMax) {
				tempMax = sum2;
				max_p = cal_indexes_range + 1;
			}
			if (sum3 > tempMax) {
				tempMax = sum3;
				max_p = cal_indexes_range + 2;
			}
			// sum of colors
			red += ((int)pix_1.red) + ((int)pix_2.red) + ((int)pix_3.red); // Instead of call sum_pixels_by_weight
			green += ((int)pix_1.green) + ((int)pix_2.green) + ((int)pix_3.green); // Instead of call sum_pixels_by_weight
			blue += ((int)pix_1.blue) + ((int)pix_2.blue) + ((int)pix_3.blue); // Instead of call sum_pixels_by_weight

			// FILTER multiple
			red += ((int)src[min_p].red) * (-1); // Instead of call sum_pixels_by_weight
			green += ((int)src[min_p].green) * (-1); // Instead of call sum_pixels_by_weight
			blue += ((int)src[min_p].blue) * (-1); // Instead of call sum_pixels_by_weight

			red += ((int)src[max_p].red) * (-1); // Instead of call sum_pixels_by_weight
			green += ((int)src[max_p].green) * (-1); // Instead of call sum_pixels_by_weight
			blue += ((int)src[max_p].blue) * (-1); // Instead of call sum_pixels_by_weight

			the_pixel.red = red / kernelScale;
			the_pixel.blue = blue / kernelScale;
			the_pixel.green = green / kernelScale;

			dst[my_index] = the_pixel;
		}
	}
}

// copy pixels from image to src.
void copyPixels(pixel* src, Image* image) {
	// The keyword register hints to compiler that a given variable can be put in a register

	register pixel* temp = (pixel*)(image->data); // temp pointer to image data and convert to pixel*
	register int multRowN, sumColRow;
	register int row, col;
	for (row = 0; row < m; row++) {
		multRowN = row * n; // calculate one time here
		for (col = 0; col < n; col++) {
			sumColRow = multRowN + col; // calculate only one time
			src[sumColRow] = temp[sumColRow]; // less calls
		}
	}
}

void doConvolution(Image* image, int kernelScale, bool filter) {

	pixel* pixelsImg = malloc(m * n * sizeof(pixel));
	copyPixels(pixelsImg,image); // pass the image and src to copy from image to src values.

	if (filter) { // if we need filter for picture
		blurFilter(kernelScale, pixelsImg,image); // blur with filter
	}
	else {
		if (kernelScale != 1) // cases: (7 or 9), level 1 - blur the picture
			blur(kernelScale, pixelsImg, image); // only parameters blur need
		else // kernal scale is 1 so level 2 is to sharp the picture.
		sharpen(kernelScale, pixelsImg, image);
	}

	free(pixelsImg);
}

void myfunction(Image* image, char* srcImgpName, char* blurRsltImgName, char* sharpRsltImgName, char* filteredBlurRsltImgName, char* filteredSharpRsltImgName, char flag) {
	if (flag == '1') {
		// blur image
		doConvolution(image, 9, false);

		// write result image to file
		writeBMP(image, srcImgpName, blurRsltImgName);

		// sharpen the resulting image
		doConvolution(image, 1, false);

		// write result image to file
		writeBMP(image, srcImgpName, sharpRsltImgName);
	}
	else {
		// apply extermum filtered kernel to blur image
		doConvolution(image, 7, true);

		// write result image to file
		writeBMP(image, srcImgpName, filteredBlurRsltImgName);

		// sharpen the resulting image
		doConvolution(image, 1, false);

		// write result image to file
		writeBMP(image, srcImgpName, filteredSharpRsltImgName);
	}
}
