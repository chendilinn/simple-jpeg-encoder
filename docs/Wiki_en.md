# Demystifying JPEG: Step-by-Step Implementation of a JPEG Encoder in C

Hey, fellow programmers and image enthusiasts! You see countless JPEG images online every day. Have you ever wondered how they become so small yet still look good? JPEG is a very popular image compression format, widely used on the web and in digital cameras.

Today, we will no longer just be users of JPEG, but its "creators"! We will explore the basic principles of JPEG encoding together, and by analyzing a simple JPEG encoder written in C (`jpeg_encoder.c`), you will experience this process firsthand. Ready? Let's embark on this wonderful journey!

## Core Steps of JPEG Encoding: The Magic of Compression

JPEG compression is not a one-step process; it's like a precise assembly line involving several key steps. Our goal is to minimize file size while maintaining the visual quality of the image as much as possible. Here are the main stages of this pipeline:

### 1. Color Space Conversion (RGB -> YCbCr): Focusing on Human Eye Sensitivity

Images we usually encounter are in RGB format (Red, Green, Blue). However, the human eye is more sensitive to **luminance** (the brightness of an image) than to **chrominance** (the color information of an image). JPEG leverages this fact.

*   **Why convert?** To separate brightness and color, allowing us to compress less sensitive color information to a greater extent.
*   **What is YCbCr?**
    *   `Y`: Luminance - represents brightness.
    *   `Cb`: Chrominance Blue - represents the difference between the blue component and luminance.
    *   `Cr`: Chrominance Red - represents the difference between the red component and luminance.

**Code Implementation (`rgb_to_ycbcr_and_subsample` function):**

In our C code, the `rgb_to_ycbcr_and_subsample` function is responsible for this conversion. It iterates through each pixel and applies the standard conversion formulas:

```c
// Iterate through each pixel for conversion
for (uint16_t y = 0; y < encoder->height; ++y) {
    for (uint16_t x = 0; x < encoder->width; ++x) {
        // ... Get R, G, B values ...

        // Use standard JPEG RGB to YCbCr conversion formulas
        double Y_f =    0.299   * R + 0.587   * G + 0.114   * B;
        double Cb_f = - 0.168736* R - 0.331264* G + 0.5     * B + 128.0; // Cb/Cr needs +128 offset
        double Cr_f =   0.5     * R - 0.418688* G - 0.081312* B + 128.0;

        // Store Y component (clamp to 0-255)
        encoder->y_data[y_idx] = clamp_byte((int)round(Y_f));

        // ... (Chroma subsampling code, see next step) ...
    }
}
// YCbCr data generated, free original RGB data memory
free(encoder->rgb_data);
encoder->rgb_data = NULL;
```
This code reads the RGB values of each pixel, calculates the corresponding Y, Cb, and Cr values using mathematical formulas, and stores the Y values. Note that 128 is added to Cb and Cr after calculation to adjust the range from potentially negative values to 0-255. After conversion, the original RGB data is freed.

### 2. Chroma Subsampling: Discarding Partial Color Information

Since the human eye is less sensitive to chrominance, can we store less chrominance information? Yes! This is chroma subsampling.

*   **Principle:** For a small area (e.g., a 2x2 pixel block), we can retain only an average chrominance value (Cb and Cr), but keep the luminance values (Y) for all pixels.
*   **Common Mode (4:2:0):** This is the most common mode. It means that for every 4 Y values (in a 2x2 block), we store only 1 Cb value and 1 Cr value. This reduces the amount of chrominance data by 75%!

**Code Implementation (`rgb_to_ycbcr_and_subsample` function):**

In the same function, immediately after calculating the Y values, the code performs subsampling for Cb and Cr:

```c
// Perform chroma subsampling for Cb and Cr (calculate and store only at the top-left pixel of a 2x2 block)
if ((y % 2 == 0) && (x % 2 == 0)) {
    size_t cbcr_idx = (y / 2) * encoder->cbcr_width + (x / 2); // Calculate index in Cb/Cr plane

    // Simple averaging for Cb/Cr values in a 2x2 region
    // ... (Calculate average Cb, Cr values cb_avg, cr_avg for the 2x2 region) ...
    double cb_sum = 0, cr_sum = 0;
    int count = 0;
    // ... (Loop through 2x2 region, accumulate Cb, Cr and count) ...

    // Store averaged Cb/Cr values (clamp to 0-255)
    encoder->cb_data[cbcr_idx] = clamp_byte((int)round(cb_sum / count));
    encoder->cr_data[cbcr_idx] = clamp_byte((int)round(cr_sum / count));
}
```
This code executes only at even rows and columns (`y % 2 == 0 && x % 2 == 0`). It calculates the average Cb and Cr values for all pixels within the current 2x2 neighborhood and then stores this average value in the reduced-size `cb_data` and `cr_data` arrays. This implements 4:2:0 subsampling.

### 3. Blocking: Divide and Conquer

JPEG does not process the entire image at once; instead, it divides it into many small 8x8 pixel blocks.

*   **Why blocking?** Subsequent DCT transformations are performed on these 8x8 blocks. This makes calculations more manageable and better suited for processing local image features.

**Code Implementation (`encode_image` function):**

Although there isn't a specific function called `blocking`, the blocking operation is reflected in how data blocks are extracted within the main encoding loop (`encode_image`). The code iterates through what are called MCUs (Minimum Coded Units, which typically contain 4 Y blocks, 1 Cb block, and 1 Cr block in 4:2:0 mode) and extracts 8x8 pixel blocks when processing each component:

```c
// Inside the MCU loop of the encode_image function...
// --- Process 4 8x8 blocks of the Y component ---
for (int v_samp = 0; v_samp < y_comp->v_samp_factor; ++v_samp) {
    for (int h_samp = 0; h_samp < y_comp->h_samp_factor; ++h_samp) {
        // 1. Extract the current 8x8 block from the Y data plane
        // ... (Calculate block start coordinates block_y_start, block_x_start) ...
        for(int y=0; y<BLOCK_SIZE; ++y) {
            for(int x=0; x<BLOCK_SIZE; ++x) {
                 // ... (Calculate source pixel coordinates src_y, src_x, handle boundaries) ...
                 input_block[y*BLOCK_SIZE + x] = encoder->y_data[src_y * encoder->width + src_x];
                 // ... (Boundary handling) ...
            }
        }
        // Next, perform DCT, Quantization, Huffman on input_block ...
    }
}
// --- Similarly process 1 8x8 block for Cb and Cr ---
// ...
```
Here, `input_block` is the 8x8 pixel block we extracted (from Y, Cb, or Cr data), and subsequent processing will be based on this small block. The `BLOCK_SIZE` constant is defined as 8.

### 4. Discrete Cosine Transform (DCT): From Pixels to Frequencies

This is the core mathematical transformation in JPEG compression.

*   **Role of DCT:** It transforms an 8x8 pixel block from the **spatial domain** (representing pixel positions and values) to the **frequency domain**. After transformation, 64 DCT coefficients are obtained.
    *   The top-left coefficient (DC coefficient) represents the **average** brightness/chrominance of the block (low-frequency information).
    *   Coefficients further to the bottom-right represent faster-changing details in the image (high-frequency information, such as edges and textures).
*   **Key Feature:** For natural images, most of the energy (important visual information) is concentrated in a few low-frequency coefficients in the top-left corner, while high-frequency coefficients in the bottom-right are usually small or close to zero. This lays the foundation for subsequent compression.
*   The one-dimensional DCT transformation formula is as follows (two-dimensional DCT is essentially performing one-dimensional DCT twice, first on each row, then on the coefficients after row DCT transformation by column):

$$ X_k = \sum_{n=0}^{N-1} x_n \cdot \cos\left(\frac{\pi}{N} \left(n + \frac{1}{2}\right) k\right), \quad k = 0, 1, 2, \ldots, N-1 $$

**Code Implementation (`perform_fdct` function):**

The `perform_fdct` function implements the forward DCT. It first performs a "Level Shift" on the input 8x8 block (subtracting 128 from values in the 0-255 range to center them around 0), and then performs the DCT calculation:

```c
// --- Forward DCT (FDCT) ---
// Performs FDCT on the input 8x8 pixel block, outputs 8x8 DCT coefficients
void perform_fdct(const uint8_t *input_block, int16_t *output_coeffs) {
    double block[BLOCK_SIZE][BLOCK_SIZE]; // Used to store level-shifted pixel values
    double temp[BLOCK_SIZE][BLOCK_SIZE];  // Intermediate calculation results

    // 1. Level Shift: Convert pixel values from 0-255 to -128 to 127 by subtracting 128
    for (int y = 0; y < BLOCK_SIZE; y++) {
        for (int x = 0; x < BLOCK_SIZE; x++) {
            block[y][x] = (double)input_block[y * BLOCK_SIZE + x] - 128.0;
        }
    }

    // 2. Perform two one-dimensional FDCTs (first on rows, then on columns)
    // 2a. Perform one-dimensional FDCT on each row
    for (int y = 0; y < BLOCK_SIZE; y++) {
        for (int u = 0; u < BLOCK_SIZE; u++) { // u is frequency index
            double sum = 0.0;
            for (int x = 0; x < BLOCK_SIZE; x++) { // x is spatial index
                sum += block[y][x] * G_FDCT_COS[x][u]; // Use pre-calculated cosine values
            }
            double c_u = (u == 0) ? 1.0 / sqrt(2.0) : 1.0; // Normalization factor C(u) for DCT-II
            temp[y][u] = c_u * sum / sqrt(BLOCK_SIZE / 2.0); // Denominator sqrt(N/2) for standard DCT-II normalization
        }
    }

    // 2b. Perform one-dimensional FDCT on each column of the intermediate result
    for (int u = 0; u < BLOCK_SIZE; u++) { // u is now column index (frequency)
        for (int v = 0; v < BLOCK_SIZE; v++) { // v is row index (frequency)
            double sum = 0.0;
            for (int y = 0; y < BLOCK_SIZE; y++) { // y is row index (spatial)
                sum += temp[y][u] * G_FDCT_COS[y][v]; // Use pre-calculated cosine values
            }
            double c_v = (v == 0) ? 1.0 / sqrt(2.0) : 1.0; // Normalization factor C(v) for DCT-II
            // Round the final DCT coefficient and store as a 16-bit integer
            output_coeffs[v * BLOCK_SIZE + u] = (int16_t)round(c_v * sum / sqrt(BLOCK_SIZE / 2.0));
        }
    }
}

// Initialize cosine values in main function:
init_fdct_cos(); // Pre-calculate FDCT cosine values

// Pre-calculate FDCT cosine values (same as IDCT)
void init_fdct_cos() {
    for (int i = 0; i < BLOCK_SIZE; i++) {
        for (int j = 0; j < BLOCK_SIZE; j++) {
            G_FDCT_COS[i][j] = cos((2.0 * i + 1.0) * j * M_PI / (2.0 * BLOCK_SIZE));
        }
    }
}
```
The code implements two-dimensional DCT by performing one-dimensional DCT twice (first on rows, then on columns). It uses pre-calculated cosine values (`G_FDCT_COS`) to speed up calculations. The `output_coeffs` array contains the 64 DCT coefficients. **Note:** The pre-calculated FDCT cosine values in the code are stored column-wise. Let's observe the data after DCT compared to the original data:

![afterDCT](https://i-blog.csdnimg.cn/direct/4968970e2ddb42bd90318bf9cc81cc81.png#pic_center =360x)

As you can see, after the DCT transformation, most of the important image information is concentrated in the top-left corner. At this point, there is no compression yet, because the original signal can still be fully restored from the DCT data. The next step, quantization, is key to compression efficiency.

### 5. Quantization: The Main Lossy Compression Step

This is the crucial step in JPEG's **lossy compression**, and it determines the final image quality and file size.

*   **Principle:** We use a "Quantization Table" to reduce the precision of the DCT coefficients. This table is also 8x8, and the larger the values in the table, the more severely the corresponding DCT coefficients are "compressed" (more precision is lost). Typically, high-frequency coefficients (corresponding to the bottom-right of the table) are assigned larger quantization values because the human eye is less sensitive to the loss of high-frequency details.
*   **How it works:** Each value in the DCT coefficient matrix is divided by the corresponding value in the quantization table and then rounded to the nearest integer.
*   **Result:** Many high-frequency coefficients become 0 after quantization, and the precision of low-frequency coefficients is also reduced. This is where information is lost, but it's also the source of compression efficiency.

**Code Implementation (`quantize` function and `scale_quant_table` function):**

The `quantize` function performs this division and rounding operation:

```c
// Quantize an 8x8 block of DCT coefficients
void quantize(int16_t *coeffs, const uint16_t *quant_table) {
    for (int i = 0; i < BLOCK_SIZE * BLOCK_SIZE; ++i) {
        // Quantization: Divide DCT coefficient by corresponding value in quantization table and round
        double div = (double)coeffs[i] / quant_table[i];
        coeffs[i] = (int16_t)floor(div + 0.5); // Simple rounding method
    }
}
```
How is this `quant_table` derived? It's obtained by scaling the standard quantization tables (defined as `STD_LUM_QUANT_TBL` and `STD_CHROM_QUANT_TBL` in the code) based on the user-specified "Quality Factor". The `scale_quant_table` function is responsible for this scaling:

```c
// Scale standard quantization table based on quality factor
void scale_quant_table(const uint8_t *std_table, uint16_t *output_table, int quality_factor) {
    // ... (Calculate scale_factor based on quality_factor) ...
    int scale_factor;
    if (quality_factor < 50) {
        scale_factor = 5000 / quality_factor;
    } else {
        scale_factor = 200 - 2 * quality_factor;
    }
    for (int i = 0; i < 64; ++i) {
        int temp = ((int)std_table[i] * scale_factor + 50) / 100; // Scale and round
        // ... (Ensure quantized value is between 1 and 255) ...
        if (temp <= 0) temp = 1;
        if (temp > 255) temp = 255;
        output_table[i] = (uint16_t)temp;
    }
}

// Called in main function:
scale_quant_table(STD_LUM_QUANT_TBL, encoder->quant_tables[0], encoder->quality_factor); // Scale luminance quantization table
scale_quant_table(STD_CHROM_QUANT_TBL, encoder->quant_tables[1], encoder->quality_factor); // Scale chrominance quantization table
```
The higher the quality factor, the smaller the `scale_factor`, and the closer the quantization table values are to the standard values (or smaller), retaining more precision, resulting in better image quality but a larger file. The opposite is also true.

Let's observe the 8x8 data after DCT and quantization:
![afterquant](https://i-blog.csdnimg.cn/direct/d4e28102762541e68e5d9fe353194bbf.png#pic_center =360x)

As you can see, after quantization, most of the data becomes 0, which is the source of compression ratio.

### 6. Zigzag Scan: Preparing for Entropy Encoding

After quantization, the 8x8 coefficient matrix often contains many zeros, especially in the bottom-right corner (high-frequency region). To compress these zeros more effectively, we don't read the coefficients in a conventional row or column order, but instead use a Zigzag scan.

*   **Why Zigzag?** It rearranges the 8x8 matrix into a one-dimensional sequence of 64 elements. Since energy is concentrated in the top-left corner (low frequencies), the Zigzag scan tends to place non-zero coefficients first, followed by a long string of zeros.

Zigzag scan order is as follows:

![ZigZag](https://i-blog.csdnimg.cn/direct/40bea984f4874dbb8614f221b4b0ee97.png#pic_center)

**Code Implementation (using `G_ZIGZAG` array):**

The code uses a pre-defined `G_ZIGZAG` array to specify the scan order. When encoding AC coefficients in the `encode_image` function, it reads the quantized coefficients in this order:

```c
// In the AC coefficient encoding loop of the encode_image function...
for (int k = 1; k < 64; ++k) { // Start from the 1st AC coefficient (in Zigzag order)
    int coeff = quant_coeffs[G_ZIGZAG[k]]; // Get the next coefficient in Zigzag scan
    if (coeff == 0) {
        zero_run++; // If it's 0, increment the count of consecutive zeros
    } else {
        // ... (Process non-zero coefficient and previous zero_run) ...
    }
}
```
`G_ZIGZAG[k]` gives the index in the original 8x8 matrix (numbered row by row, then column by column) of the `k`-th coefficient to be accessed.

### 7. Entropy Encoding: The Final Lossless Compression

This is the last step in the pipeline, aiming for **lossless compression** of the coefficient sequence after Zigzag scanning, further reducing file size. JPEG primarily uses Huffman Coding.

*   **Principle:**
    *   **DC Coefficient:** For the DC coefficient of each block (the average value of the block), we don't encode its value directly, but rather the **difference** between its value and the DC coefficient of the **previous** block (Differential Coding). Since the average brightness/chrominance of adjacent blocks is usually very similar, this difference will be small and easier to encode.
    *   **AC Coefficients:** For AC coefficients (the 1st to 63rd coefficients after Zigzag scanning), we use **Run-Length Encoding (RLE)**. Instead of storing a long string of zeros directly, we use a symbol to represent "(N consecutive zeros, followed by a non-zero coefficient V)". For example, `(5, -3)` means there are 5 zeros followed by a coefficient with a value of -3.
    *   **Huffman Coding:** The size (Category/Size) of the DC difference, the RLE symbol for AC coefficients, and the value of non-zero coefficients (how many bits are needed to represent them) are all encoded using pre-defined Huffman tables. Huffman coding assigns shorter bit strings to frequently occurring symbols and longer bit strings to less frequent symbols, thereby achieving compression. JPEG has standard Huffman tables specifically defined for DC and AC coefficients (as well as luminance/chrominance).

**Code Implementation (`encode_image`, `build_huffman_code_tables`, `write_bits`, etc.):**

This logic is primarily implemented in the Huffman encoding part of the `encode_image` function:

```c
// In encode_image...

// 4a. Encode DC coefficient (Differential Coding)
int dc_diff = quant_coeffs[0] - y_comp->dc_pred; // Calculate difference
y_comp->dc_pred = quant_coeffs[0];             // Update prediction value
int dc_category = get_category(dc_diff);       // Get Category of the difference
// Write Huffman code for DC coefficient Category
write_bits(&encoder->bitstream, dc_table_y->code[dc_category], dc_table_y->length[dc_category]);
// If Category > 0, write VLI encoding of the difference
if (dc_category > 0) {
    write_bits(&encoder->bitstream, get_vl_code(dc_diff, dc_category), dc_category);
}

// 4b. Encode AC coefficients (Run-Length Encoding RLE + Huffman)
int zero_run = 0; // Count of consecutive zeros
for (int k = 1; k < 64; ++k) {
    int coeff = quant_coeffs[G_ZIGZAG[k]];
    if (coeff == 0) {
        zero_run++;
    } else {
        // Handle accumulated zeros (write ZRL 0xF0 if zero_run >= 16)
        while (zero_run >= 16) { write_bits(&encoder->bitstream, ac_table_y->code[0xF0], ac_table_y->length[0xF0]); zero_run -= 16; }
        // Get Category (Size) of the non-zero coefficient
        int ac_category = get_category(coeff);
        // Construct RLE symbol: (run << 4) | size
        uint8_t rle_symbol = (zero_run << 4) | ac_category;
        // Write Huffman code for RLE symbol
        write_bits(&encoder->bitstream, ac_table_y->code[rle_symbol], ac_table_y->length[rle_symbol]);
        // Write VLI encoding of the non-zero coefficient
        write_bits(&encoder->bitstream, get_vl_code(coeff, ac_category), ac_category);
        zero_run = 0; // Reset
    }
}
// Handle zeros at the end of the block (write EOB 0x00)
if (zero_run > 0) { write_bits(&encoder->bitstream, ac_table_y->code[0x00], ac_table_y->length[0x00]); }
```

*   `get_category` function calculates how many bits are needed to represent a value (Category/Size).
*   `get_vl_code` function obtains the Variable-Length Integer (VLI) encoding of the coefficient value.
*   `build_huffman_code_tables` (called in `main`) pre-builds the lookup table `HuffmanCodeTable` based on standard tables (e.g., `STD_DC_LUM_BITS`, `STD_AC_LUM_VAL`), storing the Huffman code and length for each symbol.
*   `write_bits` function is responsible for writing the calculated Huffman codes and VLI values, bit by bit, to the output file's bitstream. It also handles JPEG's specific `0xFF` byte stuffing (`0xFF` must be followed by `0x00`).

### 8. JPEG File Format: Organizing Data

Finally, all this processed data, along with various metadata (such as image dimensions, quantization tables used, Huffman tables, etc.), needs to be organized according to the JPEG standard format and written to the final `.jpg` file. This involves various "Markers", which start with `0xFF` followed by a code, used to identify different types of data segments.

**Code Implementation (in `main` function: `write_marker`, `write_APP0`, `write_DQT`, `write_SOF0`, `write_DHT`, `write_SOS`):**

Our `main` function writes a series of marker segments before calling `encode_image` to process pixel data:

```c
// In main function...
write_marker(0xD8, encoder->bitstream.fp); // SOI (Start of Image)
write_APP0(encoder->bitstream.fp);         // APP0 (JFIF Header - contains version, density, etc.)
write_DQT(&encoder);                      // DQT (Define Quantization Table - writes our scaled quantization tables)
write_SOF0(&encoder);                     // SOF0 (Start of Frame - writes image dimensions, color components, sampling factors, etc.)
write_DHT(&encoder);                      // DHT (Define Huffman Table - writes used Huffman tables)
write_SOS(&encoder);                      // SOS (Start of Scan - scan starts, followed immediately by compressed image data)

// ... Call encode_image(&encoder) to write compressed data ...

write_marker(0xD9, encoder->bitstream.fp); // EOI (End of Image)
```
These `write_` functions use helper functions like `write_byte`, `write_word`, `write_marker` to write information to the file according to the JPEG specification.

### 9. Compression Effect Comparison

Let's compare the image quality and size changes with different quality factors. The original image is as follows (2.4MB):

![Original Image](https://i-blog.csdnimg.cn/direct/4f047e7b0b5248a09dce3e603f65de61.png#pic_center =560x)

Quality factor = 10, output size is 36KB. The compressed image is as follows (you can see obvious blurring):

![Quality Factor = 10](https://i-blog.csdnimg.cn/direct/7c96548bcead4f81bda22646bcd5db15.jpeg#pic_center =560x)

Quality factor = 50, output size is 109KB. The compressed image is as follows (you can see a significant improvement in image quality compared to quality factor 10, with smaller quantization coefficients and higher fidelity):

![Quality Factor = 50](Quality_50.jpeg)

### 10. How to Compile and Use This Encoder

The download links for the code and test images are available on the homepage. You can download and compile to test.

**Compilation Instructions:**

Assuming you have saved the code as `jpeg_encoder.c` and you have a C compiler (like GCC). Open your terminal or command line, navigate to the directory where the code is located, and then run:

```bash
gcc jpeg_encoder.c -o jpeg_encoder -lm
```
*   `gcc`: Invokes the GCC compiler.
*   `jpeg_encoder.c`: Your source code file.
*   `-o jpeg_encoder`: Specifies the output executable file name as `jpeg_encoder`.
*   `-lm`: **Very important!** Because the code uses functions from `math.h` (like `cos`, `round`, `floor`, `sqrt`), you need to link the math library.

If all goes well, an executable file named `jpeg_encoder` will appear in your current directory.

**Usage:**

This encoder requires three arguments:

1.  **Input File:** A PPM (Portable Pixmap) format image file. PPM is a simple, uncompressed image format, well-suited as input for our encoder. You need to prepare a `.ppm` file first (ensure it's P6 binary format, 8-bit depth). Many online tools can convert other formats to PPM.
2.  **Output File:** The name of the JPEG file you want to generate, e.g., `output.jpg`.
3.  **Quality Factor:** An integer between 1 and 100. A higher number means better image quality and a larger file. Common values are between 75 and 95.

Run the command as follows:

```bash
./jpeg_encoder input.ppm output.jpg 85
```
This command will read `input.ppm`, encode it with a quality factor of 85, and save the result as `output.jpg`, which can be viewed on Windows/Linux.

**Usage Notes/Considerations:**

*   **Input Format:** This encoder **only accepts PPM P6 format** (8-bit binary RGB) input.
*   **Sampling Mode:** It fixedly uses **4:2:0** chroma subsampling.
*   **Baseline DCT:** Implements the most basic JPEG (Baseline DCT).
*   **Standard Tables:** Uses quantization tables and Huffman tables defined in the JPEG standard.
*   **Error Handling:** The code includes some basic error checking (e.g., file opening failure), but may not be robust enough.

## Conclusion

We have completed an in-depth exploration into the inner workings of JPEG encoding. We've learned about color conversion from RGB to YCbCr, chroma subsampling that leverages human eye characteristics, 8x8 blocking for manageable processing, the magical DCT frequency transformation, quantization that determines quality and size, the clever Zigzag scan, and finally, Huffman entropy encoding. What's even better is that we've seen how these steps correspond one-to-one in the actual C code (`jpeg_encoder.c`).

Although this encoder is simple, it covers the most core principles of JPEG.
