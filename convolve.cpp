#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>
#include <cmath>
#include <string>

/*  Standard sample rate in Hz  */
#define SAMPLE_RATE       44100.0

/*  Standard sample size in bits  */
#define BITS_PER_SAMPLE   16

/*  Standard sample size in bytes  */		
#define BYTES_PER_SAMPLE  (BITS_PER_SAMPLE/8)

/*  Number of channels  */
#define MONOPHONIC        1
#define STEREOPHONIC      2

//struct WavHeader {
//    //from http://soundfile.sapp.org/doc/WaveFormat/
//    //RIFF Chunk Descriptor
//    char chunkID[4];
//    uint32_t chunkSize;
//    char format[4];

//    //fmt sub-chunk
//    char subchunk1ID[4];
//    uint32_t subchunk1Size;
//    uint16_t audioFormat;
//    uint16_t numChannels;
//    uint32_t sampleRate;
//    uint32_t byteRate;
//    uint16_t blockAlign;
//    uint16_t bitsPerSample;

//    //data sub-chunk
//    char subchunk2ID[4];
//    uint32_t subchunk2Size;
//};

//NOTE: according to http://soundfile.sapp.org/doc/WaveFormat/ " 16-bit samples are stored as 2's-complement signed integers, ranging from -32768 to 32767. "

//convert a 16-bit 2's-complement signed integer to a float scaled to -1.0 to +1.0
float sampleToFloat(int16_t value) {
   int16_t maxVal = 32767;

   float result = value / maxVal;

   return result;
}

//convert a float scaled to -1.0 to +1.0 to a 16-bit 2's-complement signed integer
int16_t floatToSample(float value) {
   int16_t maxVal = 32767;

   //this is how testtone.c does it
   int16_t result = rint(value * maxVal);

   return result;
}

//convert the vector of wav samples to floats
std::vector<float> samplesToFloats(std::vector<int16_t> samples) {
   std::vector<float> floatSamples;
   for (int16_t sample : samples) {
       float convertedSample = sampleToFloat(sample);
       floatSamples.push_back(convertedSample);
   }
   return floatSamples;
}

//convert the vector of floats into wav samples
std::vector<int16_t> floatsToSamples(std::vector<float> floatSamples) {
   std::vector<int16_t> samples;
   for (float sample : floatSamples) {
       int16_t convertedSample = floatToSample(sample);
       samples.push_back(convertedSample);
   }
   return samples;
}



std::vector<int16_t> readWavFile(char *filename) {
    // FILE *inputFileStream = fopen(filename, "rb");
    // if (inputFileStream == NULL) {
    //     fprintf(stdout, "File %s cannot be opened for reading\n", filename);
    //     return std::vector<int16_t>();
    // }

    // WavHeader header;

    // //read header
    // size_t bytesRead = fread(&header, sizeof(header), 1, inputFileStream);
    // if (bytesRead != 1) {
    //     fprintf(stdout, "Error reading WAV header.\n");
    //     fclose(inputFileStream);
    //     return std::vector<int16_t>();
    // }

    // std::cout << "subchunk2Size: " << header.subchunk2Size << std::endl; //this is way off for some reason??

    // //store samples
    // int16_t* audioSamplesArray = new int16_t[header.subchunk2Size / sizeof(int16_t)];
    // bytesRead = fread(audioSamplesArray, sizeof(int16_t), header.subchunk2Size / sizeof(int16_t), inputFileStream);
    // if (bytesRead != header.subchunk2Size / sizeof(int16_t)) {
    //     fprintf(stdout, "Error reading audio sample data.\n");
    //     delete[] audioSamplesArray;
    //     fclose(inputFileStream);
    //     return std::vector<int16_t>();
    // }

    // //sample vector
    // std::vector<int16_t> audioSamples;
    // for (size_t i = 0; i < header.subchunk2Size / sizeof(int16_t); ++i) {
    //     audioSamples.push_back(audioSamplesArray[i]);
    // }

    // //Check if the file is mono, 16-bit, 44.1 kHz
    // std::cout << header.numChannels << " " << header.bitsPerSample << " " << header.sampleRate << " " << header.subchunk2Size << std::endl;
    // if (header.numChannels != 1 || header.bitsPerSample != 16 || header.sampleRate != 44100) {
    //     std::cout << "Unsupported WAV file format. Must be mono, 16-bit, 44.1 kHz." << std::endl;
    //     return std::vector<int16_t>();
    // }
    
    // fclose(inputFileStream);
    // return audioSamples;

    // this isn't working, and apparently we're supposed to use fstream in c++ anyway https://stackoverflow.com/questions/25903590/c-vs-c-file-handling

    std::ifstream file(filename, std::ios::binary);

    if (!file.is_open()) {
        fprintf(stdout, "File %s cannot be opened for reading\n", filename);
        return std::vector<int16_t>();
    }

    //read header
    char header[44];
    file.read(header, 44);

    //Check if the file is mono, 16-bit, 44.1 kHz
    //header[22] is numchannels, 34 is bits per sample, 24 and 26 is sample rate
    if (header[22] != 1 || header[34] != 16 || *reinterpret_cast<int*>(header + 24) != 44100) {
        std::cerr << "Unsupported WAV file format. Must be mono, 16-bit, 44.1 kHz." << std::endl;
        return std::vector<int16_t>();
    }

    //store samples
    std::vector<int16_t> samples;
    int16_t sample;
    while (file.read(reinterpret_cast<char*>(&sample), sizeof(int16_t))) {
        samples.push_back(sample);
    }

    file.close();
    return samples;
}


//Note: The following are from testtone.c given on D2L
size_t fwriteIntLSB(int data, FILE *stream)
{
    unsigned char array[4];

    array[3] = (unsigned char)((data >> 24) & 0xFF);
    array[2] = (unsigned char)((data >> 16) & 0xFF);
    array[1] = (unsigned char)((data >> 8) & 0xFF);
    array[0] = (unsigned char)(data & 0xFF);
    return fwrite(array, sizeof(unsigned char), 4, stream);
}

size_t fwriteShortLSB(int16_t data, FILE *stream)
{
    unsigned char array[2];

    array[1] = (unsigned char)((data >> 8) & 0xFF);
    array[0] = (unsigned char)(data & 0xFF);
    return fwrite(array, sizeof(unsigned char), 2, stream);
}

void writeWaveFileHeader(int channels, int numberSamples,
                         double outputRate, FILE *outputFile)
{
    /*  Calculate the total number of bytes for the data chunk  */
    int dataChunkSize = channels * numberSamples * BYTES_PER_SAMPLE;
	
    /*  Calculate the total number of bytes for the form size  */
    int formSize = 36 + dataChunkSize;
	
    /*  Calculate the total number of bytes per frame  */
    int16_t frameSize = channels * BYTES_PER_SAMPLE;
	
    /*  Calculate the byte rate  */
    int bytesPerSecond = (int)ceil(outputRate * frameSize);

    /*  Write header to file  */
    /*  Form container identifier  */
    fputs("RIFF", outputFile);
      
    /*  Form size  */
    fwriteIntLSB(formSize, outputFile);
      
    /*  Form container type  */
    fputs("WAVE", outputFile);

    /*  Format chunk identifier (Note: space after 't' needed)  */
    fputs("fmt ", outputFile);
      
    /*  Format chunk size (fixed at 16 bytes)  */
    fwriteIntLSB(16, outputFile);

    /*  Compression code:  1 = PCM  */
    fwriteShortLSB(1, outputFile);

    /*  Number of channels  */
    fwriteShortLSB((int16_t)channels, outputFile);

    /*  Output Sample Rate  */
    fwriteIntLSB((int)outputRate, outputFile);

    /*  Bytes per second  */
    fwriteIntLSB(bytesPerSecond, outputFile);

    /*  Block alignment (frame size)  */
    fwriteShortLSB(frameSize, outputFile);

    /*  Bits per sample  */
    fwriteShortLSB(BITS_PER_SAMPLE, outputFile);

    /*  Sound Data chunk identifier  */
    fputs("data", outputFile);

    /*  Chunk size  */
    fwriteIntLSB(dataChunkSize, outputFile);
}

//Function to write audio data to a new WAV file
void writeWavFile(char *filename, std::vector<int16_t> samples, int numberOfChannels) {
    FILE *outputFileStream = fopen(filename, "wb");
    if (outputFileStream == NULL) {
        fprintf(stdout, "File %s cannot be opened for writing\n", filename);
        return;
    }

    int numberOfSamples = samples.size();
    std::cout << "sample num: " << numberOfSamples << std::endl;

    /*  Write the WAVE file header  */
    writeWaveFileHeader(numberOfChannels, numberOfSamples,
                        SAMPLE_RATE, outputFileStream);

    //Write audio data (samples) to the file
    for (int i = 0; i < numberOfSamples; i++) {
        fwriteShortLSB(samples[i], outputFileStream);
    }

    fclose(outputFileStream);
}


int main(int argc, char *argv[]) {
    char *inputFile = NULL;
    char *IRFile = NULL;
    char *outputFile = NULL;

    if (argc != 4) {
        std::cout << "Usage: " << argv[0] << " inputfile IRfile outputfile" << std::endl;
        return 1;
    }

    inputFile = argv[1];
    IRFile = argv[2];
    outputFile = argv[3];

    //get the data section of the input file as a vector
    std::vector<int16_t> wavSamples = readWavFile(inputFile);
    if (wavSamples.empty()) {
        return 1;
    }
    //ints to floats
    //std::vector<float> floatSamples = samplesToFloats(wavSamples);

    //convolution stuff

    //floats back to ints
    //std::vector<int16_t> outputSamples = floatsToSamples(floatSamples);

    //Write the converted samples to a new WAV file
    writeWavFile(outputFile, wavSamples, 1);

    std::cout << "Reverbification complete. Output file: " << outputFile << std::endl;

    return 0;
}