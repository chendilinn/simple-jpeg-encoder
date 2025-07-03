[English](README.md) | [简体中文](README.zh-CN.md)

# A Simple JPEG Encoder in C

This is a simple JPEG encoder written in C for educational purposes. It takes a PPM (P6) format image file as input and generates a baseline DCT-based JPEG file according to a specified quality factor.

The goal of this project is to clearly demonstrate the core steps of the JPEG compression standard, with detailed comments in the code to facilitate understanding.

## Key Features

- **PPM Input**: Supports reading binary PPM (P6) format images.
- **Color Space Conversion**: Implements conversion from the RGB to YCbCr color space.
- **Chroma Subsampling**: Supports 4:2:0 chroma subsampling to reduce the data volume of color information.
- **Forward DCT (FDCT)**: Performs a Forward Discrete Cosine Transform on 8x8 pixel blocks.
- **Quantization**: Quantizes the DCT coefficients based on a user-defined quality factor (1-100).
- **Huffman Encoding**: Encodes the quantized coefficients using standard Huffman tables.
- **Standard JPEG Output**: Generates a `.jpg` file with standard JPEG markers (SOI, APP0, DQT, SOF0, DHT, SOS, EOI).

## How to Compile

You will need a C compiler, such as `gcc`. Since the code uses the math library (`math.h`), you need to link it during compilation.

```bash
gcc jpeg_encoder.c -o jpeg_encoder -lm
```

## How to Use

After compiling, you can run the encoder from the command line.

**Usage:**
```bash
./jpeg_encoder <input.ppm> <output.jpg> <quality>
```

- `<input.ppm>`: The input image file in PPM P6 format.
- `<output.jpg>`: The desired filename for the output JPEG.
- `<quality>`: An integer from 1 to 100, where 1 is the lowest quality and 100 is the highest.

**Example:**
```bash
# Encode "input.ppm" to "output.jpg" with a quality of 80
./jpeg_encoder input.ppm output.jpg 80
```

## Core JPEG Encoding Pipeline

This encoder follows the main processing steps of the baseline JPEG standard:

1.  **Read PPM File**: Parses the header and RGB pixel data from a PPM P6 file.
2.  **Color Conversion & Subsampling**:
    - Converts the RGB value of each pixel to a YCbCr value. Y represents luminance (brightness), while Cb and Cr represent chrominance (color).
    - Applies 4:2:0 chroma subsampling to the Cb and Cr components. This means for every 2x2 block of pixels, only one Cb and one Cr value are stored, effectively compressing the data.
3.  **Blocking & Level Shifting**: Divides each component (Y, Cb, Cr) of the image into 8x8 blocks and shifts the pixel values from the [0, 255] range to [-128, 127].
4.  **Forward Discrete Cosine Transform (FDCT)**: Applies the FDCT to each 8x8 block, transforming it from the spatial domain to the frequency domain. After the transformation, most of the energy is concentrated in the top-left low-frequency coefficients.
5.  **Quantization**:
    - Adjusts the standard quantization tables based on the user-specified quality factor. Lower quality results in larger quantization values.
    - Divides each value in the 8x8 FDCT coefficient matrix by the corresponding value in the quantization table and rounds the result. This is the key lossy step in JPEG compression, where many high-frequency coefficients become zero.
6.  **Zig-Zag Scan & Huffman Encoding**:
    - Scans the quantized 8x8 coefficient matrix in a zig-zag pattern to convert the 2D matrix into a 1D sequence. This groups the many zero values together, facilitating the subsequent run-length encoding.
    - **DC Coefficient Encoding**: The DC coefficient (the first one) of each block is encoded differentially (as the difference from the previous block's DC value), and the difference is then Huffman encoded.
    - **AC Coefficient Encoding**: The AC coefficients (the other 63) are run-length encoded (RLE), and the RLE symbols are then Huffman encoded.
7.  **Write JPEG File**:
    - Writes the various marker segments (like DQT for quantization tables, DHT for Huffman tables, etc.) in the standard JPEG format.
    - Writes the compressed image data to the file.
    - Ends with the EOI (End of Image) marker.

## License

This project is licensed under the [MIT License](LICENSE).
