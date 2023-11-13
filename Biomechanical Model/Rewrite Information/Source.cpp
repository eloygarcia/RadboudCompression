#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkFlipImageFilter.h"

typedef itk::Image<unsigned char, 3> ImageType;
typedef itk::ImageFileReader< ImageType > ReaderType;
typedef itk::ImageFileWriter< ImageType > WriterType;

typedef itk::FlipImageFilter< ImageType > FlipType;

#include "itkBinaryThresholdImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkAffineTransform.h"
#include "itkResampleImageFilter.h"

typedef itk::BinaryThresholdImageFilter<ImageType, ImageType> ThresholdFilterType;
typedef itk::NearestNeighborInterpolateImageFunction<ImageType> NNInterpolationType;
typedef itk::AffineTransform< double, 3 >  TransformType;
typedef itk::ResampleImageFilter<ImageType,ImageType> ResamplingType;

void Usage(char* argv[])
{
	std::cout<< " " << std::endl;
	std::cout<< "Change information (from breast.tiff to breast.nrrd + breastMask.nrrd) " << std::endl;
	std::cout<< argv[0] <<"  ImagePath OutputPath VoxelSizeX VoxelSizeY VoxelSizeZ" << std::endl;
}


int main( int argc, char* argv[])
{
	std::cout << std::endl;
	// ===================== Usage =========================
	if(argc!=6)  //or argc!=5 ???
	{
		Usage(argv);
		return EXIT_FAILURE;
	}
	//std::string dirname = argv[1];
	std::string input = argv[1];
	std::string output = argv[2];

	//std::cout << "Directorio: " << dirname << std::endl;
	std::cout << "Input: " << input << std::endl;
	std::cout << "Output: " << output << std::endl;

	ReaderType::Pointer reader = ReaderType::New();
		reader->SetFileName( input );
		try{
			reader->Update();
		} catch( itk::ExceptionObject & excp ) {
			std::cout << "Image File Reader Exception Object! " << std::endl;
			std::cout << excp << std::endl;
			return EXIT_FAILURE;
		}

	std::cout << "Image loaded ! " << std::endl;

	//itk::FixedArray< bool, 3 > flipAxes;
	//	flipAxes[0] = 1;
	//	flipAxes[1] = 1;
	//	flipAxes[2] = 0;

	//FlipType::Pointer flip = FlipType::New();
	//	flip->SetInput( reader->GetOutput());
	//	flip->SetFlipAxes( flipAxes );
	//	flip->Update();

	// ImageType::Pointer  image = flip->GetOutput();
	ImageType::Pointer image = reader->GetOutput();
	ImageType::SizeType size = reader->GetOutput()->GetLargestPossibleRegion().GetSize();
	ImageType::PointType  origin; 
			origin[0] = 0.0;
			origin[1] = 0.0;
			origin[2] = 0.0;
		image->SetOrigin( origin );
	ImageType::DirectionType direction;
			direction.SetIdentity();
		image->SetDirection( direction );
	ImageType::SpacingType spacing;
			spacing[0] = atof(argv[3]);
			spacing[1] = atof(argv[4]);
			spacing[2] = atof(argv[5]);
		image->SetSpacing( spacing );

	std::cout << "Vamos a escribir!" << std::endl;

	WriterType::Pointer writer = WriterType::New();
		writer->SetFileName( output );
		writer->SetInput( image );
		writer->UseCompressionOn();
		try{
			writer->Update();
		} catch( itk::ExceptionObject & excp ) {
			std::cout << "Writer Exception Object!" << std::endl;
			std::cout << excp << std::endl;
			return EXIT_FAILURE;
		}

	std::cout << "Fin de la escritura" << std::endl;

	/* RESAMPLING IMAGE */
	ImageType::SpacingType spacingOutput;
		spacingOutput[0] = 0.273;
		spacingOutput[1] = 0.273;
		spacingOutput[2] = 0.273;
	ImageType::SizeType sizeOutput;
		sizeOutput[0] = int( (size[0]*spacing[0]) / spacingOutput[0]);
		sizeOutput[1] = int( (size[1]*spacing[1]) / spacingOutput[1]);
		sizeOutput[2] = int( (size[2]*spacing[2]) / spacingOutput[2]);
	/*
	std::cout << "Entra Threshold" << std::endl;
	ThresholdFilterType::Pointer threshold = ThresholdFilterType::New();
		threshold->SetInput( image );
		
		threshold->SetInsideValue( 2 );
		threshold->SetOutsideValue( 1 );

		threshold->SetLowerThreshold( 1 );
		threshold->SetUpperThreshold( 3 );

		threshold->Update();

	std::cout << "Sale Threshold" << std::endl;
	*/
	std::cout << "entra resampling!" << std::endl;
	TransformType::Pointer transform = TransformType::New();
	NNInterpolationType::Pointer interpolator = NNInterpolationType::New();
	ResamplingType::Pointer resampleFilter = ResamplingType::New();
		resampleFilter->SetInput( image );
		resampleFilter->SetDefaultPixelValue( 0 );
		
		resampleFilter->SetTransform( transform );
		resampleFilter->SetInterpolator( interpolator );
		
		resampleFilter->SetSize( sizeOutput );
		resampleFilter->SetOutputSpacing( spacingOutput );
		resampleFilter->SetOutputDirection( direction );
		resampleFilter->SetOutputOrigin( origin );
		
		resampleFilter->Update();

	std::string outputBM = output.substr(0, output.length()-4);
		outputBM = outputBM + "-Resampled.nrrd";
	std::cout << outputBM << std::endl;
	
	WriterType::Pointer writeBM = WriterType::New();
		writeBM->SetFileName( outputBM );
		writeBM->SetInput( resampleFilter->GetOutput() );
		writeBM->UseCompressionOn();
		try{
			writeBM->Update();
		} catch( itk::ExceptionObject & excp ) {
			std::cout << "Writer Exception Object!" << std::endl;
			std::cout << excp << std::endl;
			return EXIT_FAILURE;
		}


	std::cout << "SUCCESSFULL!!!" << std::endl;

	return EXIT_SUCCESS;
}