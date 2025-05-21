# Mini_Project-3: Generating_SinusoidalWaves
### 程式碼解釋:
#### 引入函式庫
```js
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
```
#### 設置WAV header
```js
#define PI 3.14159265359

    // Function to write the WAV header
    void writeWavHeader(FILE *fp, double fs, int bitsPerSample, size_t N) {
        int numChannels = 1; // Mono audio
        int byteRate = fs * bitsPerSample * numChannels / 8;
        int blockAlign = bitsPerSample * numChannels / 8;

        // Write RIFF header
        fwrite("RIFF", 1, 4, fp);
        int chunkSize = 36 + (bitsPerSample / 8) * N;
        fwrite(&chunkSize, sizeof(int), 1, fp);

        // Write WAVE format
        fwrite("WAVEfmt ", 1, 8, fp);
        int subchunk1Size = 16;
        fwrite(&subchunk1Size, sizeof(int), 1, fp);

        short audioFormat = 1; // PCM format
        fwrite(&audioFormat, sizeof(short), 1, fp);
        fwrite(&numChannels, sizeof(short), 1, fp);

        int sampleRate = (int)fs;
        fwrite(&sampleRate, sizeof(int), 1, fp);

        fwrite(&byteRate, sizeof(int), 1, fp);

        fwrite(&blockAlign, sizeof(short), 1, fp);
        fwrite(&bitsPerSample, sizeof(short), 1, fp);

        // Write data subchunk
        fwrite("data", 1, 4, fp);
        fwrite(&N, sizeof(int), 1, fp);
    }
```
#### 撰寫不同弦波的函式
```js
    // Function to generate sine wave (unchanged)
    void generateSineWave(double *x_double, short *x, size_t N, double fs, double f, double A) {
        double T = 1.0 / fs;
        for (size_t n = 0; n < N; n++) {
            double tmp = A * sin(2 * PI * f * n * T);
            x[n] = (short)floor(tmp + 0.5);
            x_double[n] = tmp;
        }
    }

    // Function to generate sawtooth wave
    void generateSawtoothWave(double *x_double, short *x, size_t N, double fs, double f, double A) {
        double T = 1.0 / fs;
        for (size_t n = 0; n < N; n++) {
            double tmp = A * (2.0 * (n * T * f - floor(n * T * f + 0.5)));
            x[n] = (short)floor(tmp + 0.5);
            x_double[n] = tmp;
        }
    }

    // Function to generate triangle wave
    void generateTriangleWave(double *x_double, short *x, size_t N, double fs, double f, double A) {
        double T = 1.0 / fs;
        for (size_t n = 0; n < N; n++) {
            double tmp = A * (4.0 * fabs((n * T * f - floor(n * T * f + 0.5))) - 1.0);
            x[n] = (short)floor(tmp + 0.5);
            x_double[n] = tmp;
        }
    }

    // Function to generate square wave
    void generateSquareWave(double *x_double, short *x, size_t N, double fs, double f, double A) {
        double T = 1.0 / fs;
        for (size_t n = 0; n < N; n++) {
            double t = fmod(n * T, 1.0 / f);
            double tmp = (t < 0.5 / f) ? A : -A;
            x[n] = (short)floor(tmp + 0.5);
            x_double[n] = tmp;
        }
    }
```
#### Main Function
```js
int main(int argc, char *argv[]) {
    if (argc != 9) {
        printf("Usage: %s fs m waveform_type f A T fn.wav sqnr.txt\n", argv[0]);
        return 1;
    }

    double fs = atof(argv[1]);
    int m = atoi(argv[2]);
    char *waveformType = argv[3]; // Store the waveform type as a string
    double f = atof(argv[4]);
    double A = atof(argv[5]);
    double T = atof(argv[6]);
    char fn[1024];
    char sqnr_fn[1024];
    strcpy(fn, argv[7]);
    strcpy(sqnr_fn, argv[8]);
    short *x = NULL; // Array of the waveform (2 bytes for each sample)
    double *x_double = NULL; // Array of the waveform as double
    size_t N = (size_t)(fs * T); // Length of the waveform (samples)
    size_t n; // Sample index
    FILE *fp = NULL; // File pointer to save the waveform

    x = (short *)malloc(sizeof(short) * N);
    x_double = (double *)malloc(sizeof(double) * N);

    if (strcmp(waveformType, "sine") == 0) {
        generateSineWave(x_double, x, N, fs, f, A);
    } else if (strcmp(waveformType, "sawtooth") == 0) {
        generateSawtoothWave(x_double, x, N, fs, f, A);
    } else if (strcmp(waveformType, "triangle") == 0) {
        generateTriangleWave(x_double, x, N, fs, f, A);
    } else if (strcmp(waveformType, "square") == 0) {
        generateSquareWave(x_double, x, N, fs, f, A);
    } else {
        printf("Invalid waveform type\n");
        free(x);
        free(x_double);
        return 1;
    }

    int bitsPerSample = m;

    fp = fopen(fn, "wb"); // Open a file pointer to save a binary file
    if (!fp) { // Check if the file is opened successfully
        fprintf(stderr, "Cannot save %s\n", fn);
        exit(1); // Stop and exit the program if an error occurs
    }

    // Write the WAV header based on the selected bit depth
    writeWavHeader(fp, fs, bitsPerSample, N);

    // Write the waveform data
    if (bitsPerSample == 8) {
        for (n = 0; n < N; n++) {
            unsigned char sample = (unsigned char)((x[n] + 32768) >> 8); // Adjust for 8-bit range
            fwrite(&sample, 1, 1, fp);
        }
    } else if (bitsPerSample == 16) {
        for (n = 0; n < N; n++) {
            short sample = x[n];
            fwrite(&sample, 2, 1, fp);
        }
    } else if (bitsPerSample == 32) {
        for (n = 0; n < N; n++) {
            int sample = x[n];
            fwrite(&sample, 4, 1, fp);
        }
    } else {
        fprintf(stderr, "Unsupported bit depth: %d\n", bitsPerSample);
        exit(1);
    }

    fclose(fp);

    // Calculate SQNR
    double signal_energy = 0.0;
    double quantization_noise_energy = 1e-20;

    for (n = 0; n < N; n++) {
        signal_energy += x_double[n] * x_double[n]; // Calculate signal energy
        quantization_noise_energy += (x_double[n] - x[n]) * (x_double[n] - x[n]); // Calculate quantization noise energy
    }

    double sqnr = 10 * log10(signal_energy / quantization_noise_energy); // Calculate SQNR in dB

    // Open a text file for writing the SQNR
    FILE *sqnrFile = fopen(sqnr_fn, "w");
    if (!sqnrFile) {
        fprintf(stderr, "Cannot open %s for writing\n", sqnr_fn);
        exit(1);
    }

    // Write SQNR value to the text file with 15 decimal places
    fprintf(sqnrFile, "SQNR: %.15f\n", sqnr);

    fclose(sqnrFile); // Close the text file

    free(x); // Free the allocated memory
    free(x_double);

    return 0;
}
```
