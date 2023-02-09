#include "CudaProjection.h"
#include "CudaLibrary_kernels.cu"

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

// Constructor y Destructor

CudaProjection::CudaProjection()
{
	// initial ct image
	m_initial_size = new int[3];
	m_initial_spacing = new float[3];
	m_initial_origen = new float[3];
	numberOfPixels = 0;

    // grid
	m_grid_origen = new float[3];
	m_grid_spacing = new float[3];
	m_grid_size = new int[3];
	

	// mamo simulada
	m_simulada_size = new int[3];
	m_simulada_spacing = new float[3];
	m_simulada_origen = new float[3];
	numberOfPixels_Simulada=0;

	cudaError_t cudaStatus;
}

CudaProjection::~CudaProjection()
{
	delete[] m_initial_size;
	delete[] m_initial_spacing;
	delete[] m_initial_origen;
	delete[] m_initial_imagepointer;

	delete [] m_i_points;
	delete [] m_f_points;
	delete [] m_elements;

	delete [] m_grid_origen;
	delete [] m_grid_size;
	delete [] m_grid_spacing;

	delete [] m_flags;
	delete [] m_cumsum;
	delete [] m_correspondingElements;

	delete[] m_simulada_size;
	delete[] m_simulada_spacing;
	delete[] m_simulada_origen;
}

void CudaProjection::Initialize()
{
	/*m_mri_size[0] = m_parameters->mri_size[0];
	m_mri_size[1] = m_parameters->mri_size[1];
	m_mri_size[2] = m_parameters->mri_size[2];*/
// numberOfPixels_MRI = m_mri_size[0] * m_mri_size[1] * m_mri_size[2];

	/*m_mri_spacing[0] = m_parameters->mri_spacing[0];
	m_mri_spacing[1] = m_parameters->mri_spacing[1];
	m_mri_spacing[2] = m_parameters->mri_spacing[2];*/

	/*m_mri_origen[0] = m_parameters->mri_origen[0];
	m_mri_origen[1] = m_parameters->mri_origen[1];
	m_mri_origen[2] = m_parameters->mri_origen[2];*/

//	m_mri_imagepointer = new float[ numberOfPixels_MRI ];
//		for(int i=0; i<numberOfPixels_MRI; i++) m_mri_imagepointer[i] = m_parameters->mri_imagePointer[i];

// numberOfPoints = m_parameters->numberOfPoints;
//	m_i_points = new float[ 3*numberOfPoints ];
//		for( int i=0; i<3*numberOfPoints; i++) m_i_points[i] = m_parameters->initial_points[i];
//	m_f_points = new float[ 3*numberOfPoints ];
//		for(int i=0; i<3*numberOfPoints; i++) m_f_points[i] = m_parameters->final_points[i];
// numberOfElements = m_parameters->numberOfElements;
//	m_elements = new int[4*numberOfElements];	
//	for(int i=0; i<4*numberOfElements; i++) m_elements[i] = m_parameters->elements[i];
//
//	m_grid_origen[0] = m_parameters->grid_origen[0];
//	m_grid_origen[1] = m_parameters->grid_origen[1];
//	m_grid_origen[2] = m_parameters->grid_origen[2];
//
//	m_grid_spacing[0] = m_parameters->grid_spacing[0];
//	m_grid_spacing[1] = m_parameters->grid_spacing[1];
//	m_grid_spacing[2] = m_parameters->grid_spacing[2];
//
//	m_grid_size[0] = m_parameters->grid_size[0];
//	m_grid_size[1] = m_parameters->grid_size[1];
//	m_grid_size[2] = m_parameters->grid_size[2];
// numberOfVoxelsGrid = m_grid_size[0] * m_grid_size[1] * m_grid_size[2] ; 

	//m_flags = new int[ numberOfVoxelsGrid ];
	//	for(int i=0; i<numberOfVoxelsGrid; i++) m_flags[i] = m_parameters->flags[i];
	//m_cumsum = new int[ numberOfVoxelsGrid ];
	//	for(int i=0; i<numberOfVoxelsGrid; i++) m_cumsum[i] = m_parameters->cumsum[i];
//maximumCorrespondingElements = m_cumsum[ numberOfVoxelsGrid -1];
//	m_correspondingElements = new int[ maximumCorrespondingElements];
//	for (int i=0; i< maximumCorrespondingElements; i++)	m_correspondingElements[i] = m_parameters->correspondingElements[i];
//
//	m_simulada_size[0] = m_parameters->mamo_size[0];
//	m_simulada_size[1] = m_parameters->mamo_size[1];
//	m_simulada_size[2] = m_parameters->mamo_size[2];
//numberOfPixels_Simulada = m_simulada_size[0] * m_simulada_size[1];
//
//	m_simulada_origen[0] = m_parameters->mamo_origen[0];
//	m_simulada_origen[1] = m_parameters->mamo_origen[1];
//	m_simulada_origen[2] = m_parameters->mamo_origen[2];
//
//	m_simulada_spacing[0] = m_parameters->mamo_spacing[0];
//	m_simulada_spacing[1] = m_parameters->mamo_spacing[1];
//	m_simulada_spacing[2] = m_parameters->mamo_spacing[2];
//	
//		//m_simulada_imagepointer = new float[numberOfPixels_Simulada];
//	m_simulada_imagepointer = new unsigned short[numberOfPixels_Simulada];
//
//	m_source[0] = m_parameters->source[0];
//	m_source[1] = m_parameters->source[1];
//	m_source[2] = m_parameters->source[2];
}


void CudaProjection::SetInitialImage(int* imageSize, 
                                    float* imageSpacing,
                                    float* imageOrigin, 
                                    unsigned char* imagePointer)
{
    m_initial_size = imageSize;
numberOfPixels = (long)imageSize[0]*(long)imageSize[1]*(long)imageSize[2];
    m_initial_spacing = imageSpacing;
    m_initial_origen = imageOrigin;
    m_initial_imagepointer = imagePointer;
}

void CudaProjection::SetMesh(int numPoints, float* i_points, float* f_points, int numElements, int* elements)
{
numberOfPoints = numPoints;
    m_i_points = i_points;
    m_f_points = f_points;
numberOfElements = numElements;
    m_elements = elements;
}

void CudaProjection::SetGrid(int* gridSize, float* gridSpacing, float* gridOrigin,
                            int* flags, int* cumsum, int* correspondingElements)
{
    m_grid_size = gridSize;
numberOfVoxelsGrid = gridSize[0]*gridSize[1]*gridSize[2];
    m_grid_spacing = gridSpacing;
    m_grid_origen = gridOrigin;

    m_flags = flags;
    m_cumsum = cumsum;
    m_correspondingElements = correspondingElements;
maximumCorrespondingElements = cumsum[numberOfVoxelsGrid -1];
}

void CudaProjection::SetFinalImage(int* finalSize, float* finalSpacing, float* finalOrigin, unsigned char* finalPointer)
{
    m_simulada_size = finalSize;
    m_simulada_spacing = finalSpacing;
    m_simulada_origen = finalOrigin;
numberOfPixels_Simulada = (long)finalSize[0]*(long)finalSize[1]*(long)finalSize[2];
    m_simulada_imagepointer = finalPointer;
}
// Metodos
void CudaProjection::Update()
{
	/*
	printf("\n");
	printf("Entra en cuda\n");
	printf("\n");

// Timer !
	cudaEvent_t start, stop;
	float time;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	cudaEventRecord(start,0);
	*/

// Initialize(); // === Supongo inicializado desde el principio.
	m_simulada_imagepointer = new unsigned char[numberOfPixels_Simulada];

// Allocacion de la memoria GPU !
	// MRI !
    cudaStatus = cudaMalloc((void**)&dev_initial_size, 3*sizeof(int));
	    if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA alloc initial Size!\n");   // MRI Size
	cudaStatus = cudaMalloc((void**)&dev_initial_spacing, 3*sizeof(float));
		if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA alloc initial Spacing!\n"); // MRI Spacing
	cudaStatus = cudaMalloc((void**)&dev_initial_origen, 3*sizeof(float));
		if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA alloc initial Origen!\n");  // MRI Origen

	cudaStatus = cudaMalloc((void**)&dev_initial_imagepointer, numberOfPixels*sizeof(unsigned char));
		if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA alloc initial image pointer!\n");  // MRI imagen !!

	// Mesh & Grid !!
	cudaStatus = cudaMalloc((void**)&dev_i_points, 3*numberOfPoints*sizeof(float)); 
		if (cudaStatus != cudaSuccess) fprintf( stderr, "CUDA alloc i_point !\n"); // i_points
	cudaStatus = cudaMalloc((void**)&dev_f_points, 3*numberOfPoints*sizeof(float));
		if (cudaStatus != cudaSuccess) fprintf( stderr, "CUDA alloc f_point !\n"); // f_points

	cudaStatus = cudaMalloc((void**)&dev_elements, 4*numberOfElements*sizeof(int));
		if (cudaStatus != cudaSuccess) fprintf( stderr, "CUDA alloc elements !\n"); // elements;

	cudaStatus = cudaMalloc((void**)&dev_grid_origen, 3*sizeof(float));
		if (cudaStatus != cudaSuccess) fprintf( stderr, "CUDA alloc Grid Origen !\n"); // grid_origen
	cudaStatus = cudaMalloc((void**)&dev_grid_spacing, 3*sizeof(float));
		if (cudaStatus != cudaSuccess) fprintf( stderr, "CUDA alloc Grid spacing!\n"); // grid_Spacing
	cudaStatus = cudaMalloc((void**)&dev_grid_size, 3*sizeof(int));
		if (cudaStatus != cudaSuccess) fprintf( stderr, "CUDA alloc Grid Size !\n"); // grid_Size

	cudaStatus = cudaMalloc((void**)&dev_flags, numberOfVoxelsGrid*sizeof(int));
		if (cudaStatus != cudaSuccess) fprintf( stderr, "CUDA alloc Flags !\n"); // flags
	cudaStatus = cudaMalloc((void**)&dev_cumsum, numberOfVoxelsGrid*sizeof(int));
		if (cudaStatus != cudaSuccess) fprintf( stderr, "CUDA alloc CumSum !\n"); // cumsum
	cudaStatus = cudaMalloc((void**)&dev_correspondingElements, maximumCorrespondingElements*sizeof(int));
		if (cudaStatus != cudaSuccess) fprintf( stderr, "CUDA alloc corresponding element !\n"); // corresponding Elements !

	// Imagen Simulada !
	cudaStatus = cudaMalloc((void**)&dev_simulada_size, 3*sizeof(int));
		if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA alloc simulada Size!\n");  // Simulada Size
	cudaStatus = cudaMalloc((void**)&dev_simulada_spacing, 3*sizeof(float));
		if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA alloc simulada Spacing!\n");  // Simulada Spacing
	cudaStatus = cudaMalloc((void**)&dev_simulada_origen, 3*sizeof(float));
		if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA alloc simulada Origen!\n");  // simulada Origen

	cudaStatus = cudaMalloc((void**)&dev_simulada_imagepointer, numberOfPixels_Simulada*sizeof(unsigned char)); // u
		if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA alloc simulada image pointer!\n");  // Simulada Imagen !!
	cudaStatus = cudaMemset((void*)dev_simulada_imagepointer, 0, numberOfPixels_Simulada*sizeof(unsigned char));
		if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA alloc simulada image pointer!\n");  // Inicializaci\F3n de la imagen simulada a Zeros !!

	// Source !
//	cudaStatus = cudaMalloc((void**)&dev_source, 3*sizeof(float));
//		if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA alloc Source!\n");  // source !"


// Copia a memoria device
	// mri
	cudaStatus = cudaMemcpy(dev_initial_size, (const int*) m_initial_size, 3*sizeof(int), cudaMemcpyHostToDevice);
	// cudaStatus = cudaMemcpy(dev_mri_size, (const int*) m_mri_size, 3*sizeof(int), cudaMemcpyHostToDevice);
	    if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA memcopy initial Size Hostsimuladaev!\n");
	cudaStatus = cudaMemcpy(dev_initial_spacing, (const float*) m_initial_spacing, 3*sizeof(float), cudaMemcpyHostToDevice);
	// cudaStatus = cudaMemcpy(dev_mri_spacing, this->m_mri_spacing, 3*sizeof(float), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA memcopy initial Spacing Hostsimuladaev!\n");
	cudaStatus = cudaMemcpy(dev_initial_origen, (const float*) m_initial_origen, 3*sizeof(float), cudaMemcpyHostToDevice);
	// cudaStatus = cudaMemcpy(dev_mri_origen, this->m_mri_origen, 3*sizeof(float), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA memcopy initial Origen Hostsimuladaev!\n");
	cudaStatus = cudaMemcpy(dev_initial_imagepointer, (const unsigned char*) m_initial_imagepointer, numberOfPixels*sizeof(unsigned char), cudaMemcpyHostToDevice);
	//cudaStatus = cudaMemcpy(dev_mri_imagepointer, this->m_mri_imagepointer, numberOfPixels_MRI*sizeof(float), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA memcopy initial image pointer Hostsimuladaev!\n");

	// Mesh & Grid !!
	cudaStatus = cudaMemcpy(dev_i_points, (const float*) m_i_points, 3*numberOfPoints*sizeof(float), cudaMemcpyHostToDevice); 
		if (cudaStatus != cudaSuccess) fprintf( stderr, "CUDA memcopy i_point Hostsimuladaev !\n"); // i_points
	cudaStatus = cudaMemcpy(dev_f_points, (const float*) m_f_points,  3*numberOfPoints*sizeof(float), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) fprintf( stderr, "CUDA memcopy f_point Hostsimuladaev !\n"); // f_points

	cudaStatus = cudaMemcpy(dev_elements, (const int*) m_elements, 4*numberOfElements*sizeof(int), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) fprintf( stderr, "CUDA memcopy elements Hostsimuladaev !\n"); // elements;

	cudaStatus = cudaMemcpy(dev_grid_origen, (const float*) m_grid_origen, 3*sizeof(float), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) fprintf( stderr, "CUDA memcopy Grid Origen !\n"); // grid_origen
	cudaStatus = cudaMemcpy(dev_grid_spacing, (const float*) m_grid_spacing, 3*sizeof(float), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) fprintf( stderr, "CUDA memcopy Grid spacing!\n"); // grid_Spacing
	cudaStatus = cudaMemcpy(dev_grid_size, (const int*) m_grid_size, 3*sizeof(int), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) fprintf( stderr, "CUDA memcopy Grid Size !\n"); // grid_Size

	cudaStatus = cudaMemcpy(dev_flags, (const int*) m_flags, numberOfVoxelsGrid*sizeof(int), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) fprintf( stderr, "CUDA memcpy Flags !\n"); // flags
	cudaStatus = cudaMemcpy(dev_cumsum, (const int*) m_cumsum, numberOfVoxelsGrid*sizeof(int), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) fprintf( stderr, "CUDA memcpy CumSum !\n"); // cumsum
	cudaStatus = cudaMemcpy(dev_correspondingElements, (const int*) m_correspondingElements, maximumCorrespondingElements*sizeof(int), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) fprintf( stderr, "CUDA memcpy corresponding element !\n"); // corresponding Elements !

	// mamo simulada
	cudaStatus = cudaMemcpy(dev_simulada_size, (const int*) m_simulada_size, 3*sizeof(int), cudaMemcpyHostToDevice);
	//cudaStatus = cudaMemcpy(dev_simulada_size, this->m_simulada_size, 3*sizeof(int), cudaMemcpyHostToDevice);
	    if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA memcopy simulada Size Hostsimuladaev!\n");
	cudaStatus = cudaMemcpy(dev_simulada_spacing, (const float*) m_simulada_spacing, 3*sizeof(float), cudaMemcpyHostToDevice);
	//cudaStatus = cudaMemcpy(dev_simulada_spacing, this->m_simulada_spacing, 3*sizeof(float), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA memcopy simulada Spacing Hostsimuladaev!\n");
	cudaStatus = cudaMemcpy(dev_simulada_origen, (const float*) m_simulada_origen, 3*sizeof(float), cudaMemcpyHostToDevice);
	//cudaStatus = cudaMemcpy(dev_simulada_origen, this->m_simulada_origen, 3*sizeof(float), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA memcopy simulada Origen Hostsimuladaev!\n");
//	cudaStatus = cudaMemcpy(dev_simulada_imagepointer, this->m_simulada_imagepointer, numberOfPixels_simulada*sizeof(float), cudaMemcpyHostToDevice);
//		if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA memcopy simulada image pointer Hostsimuladaev!\n");

//	cudaStatus = cudaMemcpy(dev_source, (const float*) m_source, 3*sizeof(float), cudaMemcpyHostToDevice);
//	//cudaStatus = cudaMemcpy(dev_source, this->m_source, 3*sizeof(float), cudaMemcpyHostToDevice);
//	    if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA memcopy initial Source Hostsimuladaev!\n");


		// Checkin memory !!
		printf("\n");
		printf("initial size: [ %i,%i,%i]\n", m_initial_size[0],m_initial_size[1],m_initial_size[2]);
		printf("initial Spacing : [ %f,%f,%f]\n",m_initial_spacing[0], m_initial_spacing[1], m_initial_spacing[2]);
		printf("initial Origen : [ %f,%f,%f]\n",m_initial_origen[0], m_initial_origen[1], m_initial_origen[2]);
		printf("Number of Pixels : %lu\n",numberOfPixels);
		
		printf("simulada Size : [ %d,%d,%d]\n",m_simulada_size[0], m_simulada_size[1], m_simulada_size[2]);
		printf("simulada Spacing : [ %f,%f,%f]\n",m_simulada_spacing[0], m_simulada_spacing[1], m_simulada_spacing[2]);
		printf("simulada Origen : [ %f,%f,%f]\n",m_simulada_origen[0], m_simulada_origen[1], m_simulada_origen[2]);
		printf("Number of Pixels : %lu\n",numberOfPixels_Simulada);
		
		//printf("Source : [ %f,%f,%f]\n",m_source[0], m_source[1], m_source[2]);
		printf("\n");


		// bl = (int)ceilf((float)(numberOfPixels_Simulada/512))+1;
		bl = (unsigned int) (1+(numberOfPixels_Simulada/((long)512)));
		printf("Number of blocks: %d\n", bl);
		printf("\n"); 
// Kernel de proyecci\F3n?
		printf("Entra en el kernel\n" );
		 // fill_dos <<< bl,512 >>> (dev_simulada_imagepointer);
		
	//cudaStatus = cudaSetDevice(0);
	kernel_projection <<< bl,512 >>> (dev_initial_size, dev_initial_spacing, dev_initial_origen, dev_initial_imagepointer,
										dev_i_points, dev_f_points, dev_elements,
										dev_grid_origen, dev_grid_spacing, dev_grid_size,
										dev_flags, dev_cumsum, dev_correspondingElements,
									  dev_simulada_size, dev_simulada_spacing, dev_simulada_origen, dev_simulada_imagepointer);
									  //dev_source);								   
					   
	cudaDeviceSynchronize();

	cudaError_t error = cudaGetLastError();
	if(error!=cudaSuccess)
	{
   		fprintf(stderr,"ERROR: %s\n", cudaGetErrorString(error) );
   		exit(-1);
	}
	printf("Sale del kernel\n" );		

// Copia a memoria host
/*	int temp_initialsize[3] = {0,0,0};
	cudaStatus = cudaMemcpy( temp_initialsize, (const int*) dev_initial_size, 3*sizeof(int),  cudaMemcpyDeviceToHost);
	    if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA memcopy initial Size Dev2Host!\n");
		else m_initial_size = temp_initialsize;

	float temp_initialspacing[3] = {0.0, 0.0,0.0};
	cudaStatus = cudaMemcpy( temp_initialspacing, (const float*) dev_initial_spacing, 3*sizeof(float),  cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA memcopy initial Spacing Dev2Host!\n");
		else m_initial_spacing = temp_initialspacing;

	float temp_initialorigen[3] = {0.0,0.0,0.0};
	cudaStatus = cudaMemcpy( temp_initialorigen, (const float*) dev_initial_origen, 3*sizeof(float),  cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA memcopy initial Origen Dev2Host!\n");
		else m_initial_origen = temp_initialorigen;

	float* imageinitialpointer = new float[numberOfPixels_initial];
	cudaStatus = cudaMemcpy( imageinitialpointer, (const float*) dev_initial_imagepointer, numberOfPixels_initial*sizeof(float), cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA memcopy initial image pointer Dev2Host!\n");
		else m_initial_imagepointer = imageinitialpointer;

	int temp_simuladasize[3] = {0,0,0};
	cudaStatus = cudaMemcpy( temp_simuladasize, (const int*) dev_simulada_size, 3*sizeof(int), cudaMemcpyDeviceToHost);
	    if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA memcopy simulada Size Dev2Host!\n");
		else m_simulada_size = temp_simuladasize;

	float temp_simuladaspacing[3] = {0.0,0.0,0.0};
	cudaStatus = cudaMemcpy( temp_simuladaspacing, (const float*) dev_simulada_spacing, 3*sizeof(float), cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA memcopy simulada Spacing Dev2Host!\n");
		else m_simulada_spacing = temp_simuladaspacing;

	float temp_simuladaorigen[3] = {0.0,0.0,0.0};
	cudaStatus = cudaMemcpy( temp_simuladaorigen, (const float*) dev_simulada_origen, 3*sizeof(float), cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA memcopy simulada Origen Dev2Host!\n");
		else m_simulada_origen = temp_simuladaorigen;
*/
//	printf("va a inicializar con pixels \n" );
	//float * temp_imagesimuladapointer;
	 //float * temp_imagesimuladapointer = new float[numberOfPixels_Simulada]; // El m_... no convenci\F3 quiz\E1 porque no esta ba inicializado...

	cudaStatus = cudaMemcpy(m_simulada_imagepointer, (const char*) dev_simulada_imagepointer, numberOfPixels_Simulada*sizeof(char),  cudaMemcpyDeviceToHost);
	//cudaStatus = cudaMemcpy(m_simulada_imagepointer, (const unsigned short*)dev_simulada_imagepointer, numberOfPixels_Simulada*sizeof(unsigned short),  cudaMemcpyDeviceToHost);

//		if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA memcopy simulada image pointer Dev2Host!\n");
//		else  m_parameters->simulada_imagePointer = m_simulada_imagepointer;
	
/*  AQUI VA LA RECUPERACION DE LA IMAGEN ORIGINAL !! RECUERDALO PORQUE ESTO HABRA QUE CAMBIARLO

	printf("va a inicializar con pixels \n" );
	//float * temp_imagesimuladapointer;
	 float * temp_imagesimuladapointer = new float[numberOfPixels_Simulada]; // El m_... no convenci\F3 quiz\E1 porque no esta ba inicializado...
	//cudaStatus = cudaMemcpy(temp_imagesimuladapointer, (const float*) dev_simulada_imagepointer, numberOfPixels_Simulada*sizeof(float),  cudaMemcpyDeviceToHost);
	cudaStatus = cudaMemcpy(temp_imagesimuladapointer, (const float*)dev_simulada_imagepointer, numberOfPixels_Simulada*sizeof(float),  cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA memcopy simulada image pointer Dev2Host!\n");
		else {
			m_simulada_imagepointer = temp_imagesimuladapointer;
			m_parameters->simulada_imagePointer = temp_imagesimuladapointer;
		}


*/


/*	float temp_source[3] = {0.0,0.0,0.0};
	cudaStatus = cudaMemcpy( temp_source, (const float*) dev_source, 3*sizeof(float), cudaMemcpyDeviceToHost);
	    if (cudaStatus != cudaSuccess) fprintf(stderr, "CUDA memcopy Source Dev2Host!\n");
		else m_source = temp_source;
*/

/*
	printf("\n");
	printf("initial size: [ %i,%i,%i]\n", m_initial_size[0],m_initial_size[1],m_initial_size[2]);
	printf("initial Spacing : [ %f,%f,%f]\n",m_initial_spacing[0], m_initial_spacing[1], m_initial_spacing[2]);
	printf("initial Origen : [ %f,%f,%f]\n",m_initial_origen[0], m_initial_origen[1], m_initial_origen[2]);
	printf("initial Image Pointer : [ %f,%f,%f, ...]\n",m_initial_imagepointer[0], m_initial_imagepointer[1], m_initial_imagepointer[2]);


	printf("simulada Size : [ %d,%d,%d]\n",m_simulada_size[0], m_simulada_size[1], m_simulada_size[2]);
	printf("simulada Spacing : [ %f,%f,%f]\n",m_simulada_spacing[0], m_simulada_spacing[1], m_simulada_spacing[2]);
	printf("simulada Origen : [ %f,%f,%f]\n",m_simulada_origen[0], m_simulada_origen[1], m_simulada_origen[2]);
	printf("simulada Image Pointer : [ %f,%f,%f, ...]\n",m_simulada_imagepointer[0], m_simulada_imagepointer[1], m_simulada_imagepointer[2]);

	printf("Source : [ %f,%f,%f]\n",m_source[0], m_source[1], m_source[2]);
	printf("\n");
*/

// Liberando memoria !!
	cudaFree( (void*) dev_initial_size);				//cudaFree( temp_initialsize);
	cudaFree( (void*) dev_initial_spacing);			//cudaFree( temp_initialspacing);
	cudaFree( (void*) dev_initial_origen );			//cudaFree( temp_initialorigen);
	cudaFree( (void*) dev_initial_imagepointer );	//cudaFree( imageinitialpointer);

	cudaFree( (void*) dev_i_points );
	cudaFree( (void*) dev_f_points );
	cudaFree( (void*) dev_elements );

	cudaFree( (void*) dev_grid_origen );
	cudaFree( (void*) dev_grid_spacing );
	cudaFree( (void*) dev_grid_size );

	cudaFree( (void*) dev_flags );
	cudaFree( (void*) dev_cumsum );
	cudaFree( (void*) dev_correspondingElements );

	cudaFree( (void*) dev_simulada_size);				//cudaFree( temp_simuladasize);
	cudaFree( (void*) dev_simulada_spacing);			//cudaFree( temp_simuladaspacing);
	cudaFree( (void*) dev_simulada_origen);			//cudaFree( temp_simuladaorigen);
	cudaFree( (void*) dev_simulada_imagepointer);		//cudaFree( (void*) temp_imagesimuladapointer);

	//cudaFree( numberOfPixels_initial);		cudaFree( numberOfPixels_simulada);

//	cudaFree( (void*) dev_source);				//cudaFree( temp_source);
	cudaDeviceReset();

	// cudaFree( kernel_projection );

// Time !!
/*	cudaEventRecord(stop,0);
	//cudaEventSynchronize( stop);
	cudaEventElapsedTime( &time, start, stop);
	printf( "Time: %f ms.\n", time);

	printf("\n");
	printf("Sale de cuda\n");
	printf("\n");
	*/
}
