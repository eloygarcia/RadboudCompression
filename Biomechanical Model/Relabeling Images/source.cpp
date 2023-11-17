#include "itkImage.h"
#include "itkImageFileReader.h"
// #include "itkBinaryThresholdImageFilter.h"
#include "itkImageFileWriter.h"

typedef itk::Image<unsigned short,3> ImageType;
typedef itk::ImageFileReader<ImageType> ReaderType;
// typedef itk::BinaryThresholdImageFilter<ImageType, ImageType> ThresholdType;
typedef itk::ImageFileWriter< ImageType > WriterType;

//#include "itkEquivalenceTable.h"
#include "itkChangeLabelImageFilter.h"
typedef itk::ChangeLabelImageFilter<ImageType, ImageType> ChangingType;


int main(int argc, char* argv[]){

    std::string inputfilename = argv[1];
    std::string outputfilename = argv[2];

    // Itk Image
    ReaderType::Pointer reader = ReaderType::New();
        reader->SetFileName( inputfilename );
        try{
            reader->Update();
            std::cout << "Image Reader ok!"<< std::endl;
        } catch( itk::ExceptionObject & excp){
            std::cout << "Reader Exception Object!" << std::endl;
            std::cout << excp << std::endl;
            return EXIT_FAILURE;
        }

    ImageType::Pointer img = reader->GetOutput();

    //////
    // 1. fatty
    // 2. glandular
    // 3. skin
    // 4. muscle
    // 5. cancer
    // 6 .paddle
    std::map<short unsigned int,short unsigned int> mapping{
         {1, 1}, // Fatty Tissue
         {2, 3}, // Skin
        {29, 2}, // glandular tissue
        {33, 3}, // nipple as skin tissue
        {40, 4}, //muscle
        {88, 1}, // ligament as glandular tissue
        {95, 1}, // TDLU as glandular tissue
       {125, 1}, // duct as glandular tissue
       {150, 1}, // artery as fatty tissue
       {225, 1}, // vein as fatty tissue
       {200, 5}, // cancerous mass
       {250, 2}, // calcification as glandular tissue
        {50, 6}, // compression paddel
    };
    ChangingType::ChangeMapType map = mapping;

    ChangingType::Pointer change = ChangingType::New();
        change->SetInput( img );
        //pix1 = 40;
        //pix2 = 0;
        //change->SetChange( pix1, pix2);
        change->SetChangeMap( map );
        change->Update();

    ////
    WriterType::Pointer writer = WriterType::New();
        writer->SetFileName( outputfilename );
        writer->SetInput( change->GetOutput() );
        writer->UseCompressionOn();
        try{
            writer->Update();
            std::cout << "Image Writer OK!" << std::endl;
        } catch( itk::ExceptionObject & excp ) {
            std::cout << "Writer Exception Object!" << std::endl;
            std::cout << excp << std::endl;
            return EXIT_FAILURE;
        };

    return 0;
}