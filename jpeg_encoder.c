#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h> // 用于固定宽度的整数类型，如 uint8_t, uint16_t

// --- 常量定义 ---
#define BLOCK_SIZE 8        // JPEG 处理的基本块大小为 8x8 像素
#define MAX_COMPONENTS 3    // 最大颜色分量数 (Y, Cb, Cr)

// --- 标准表 ---

// 标准亮度量化表 (Standard Luminance Quantization Table)
const uint8_t STD_LUM_QUANT_TBL[BLOCK_SIZE * BLOCK_SIZE] = {
    16, 11, 10, 16, 24, 40, 51, 61,
    12, 12, 14, 19, 26, 58, 60, 55,
    14, 13, 16, 24, 40, 57, 69, 56,
    14, 17, 22, 29, 51, 87, 80, 62,
    18, 22, 37, 56, 68, 109, 103, 77,
    24, 35, 55, 64, 81, 104, 113, 92,
    49, 64, 78, 87, 103, 121, 120, 101,
    72, 92, 95, 98, 112, 100, 103, 99
};

// 标准色度量化表 (Standard Chrominance Quantization Table)
const uint8_t STD_CHROM_QUANT_TBL[BLOCK_SIZE * BLOCK_SIZE] = {
    17, 18, 24, 47, 99, 99, 99, 99,
    18, 21, 26, 66, 99, 99, 99, 99,
    24, 26, 56, 99, 99, 99, 99, 99,
    47, 66, 99, 99, 99, 99, 99, 99,
    99, 99, 99, 99, 99, 99, 99, 99,
    99, 99, 99, 99, 99, 99, 99, 99,
    99, 99, 99, 99, 99, 99, 99, 99,
    99, 99, 99, 99, 99, 99, 99, 99
};

// 标准 DC 亮度霍夫曼表 (Bits - 各长度码字数量, Values - 符号值)
const uint8_t STD_DC_LUM_BITS[17] = { 0, 0, 1, 5, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0 }; // bits[0] 未使用
const uint8_t STD_DC_LUM_VAL[12] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 }; // 对应上面 bits 定义的符号

// 标准 AC 亮度霍夫曼表 (Bits, Values)
const uint8_t STD_AC_LUM_BITS[17] = { 0, 0, 2, 1, 3, 3, 2, 4, 3, 5, 5, 4, 4, 0, 0, 1, 125 };
const uint8_t STD_AC_LUM_VAL[162] = { /* ... (省略具体值，与英文版相同) ... */
    0x01, 0x02, 0x03, 0x00, 0x04, 0x11, 0x05, 0x12, 0x21, 0x31, 0x41, 0x06, 0x13, 0x51, 0x61, 0x07,
    0x22, 0x71, 0x14, 0x32, 0x81, 0x91, 0xA1, 0x08, 0x23, 0x42, 0xB1, 0xC1, 0x15, 0x52, 0xD1, 0xF0,
    0x24, 0x33, 0x62, 0x72, 0x82, 0x09, 0x0A, 0x16, 0x17, 0x18, 0x19, 0x1A, 0x25, 0x26, 0x27, 0x28,
    0x29, 0x2A, 0x34, 0x35, 0x36, 0x37, 0x38, 0x39, 0x3A, 0x43, 0x44, 0x45, 0x46, 0x47, 0x48, 0x49,
    0x4A, 0x53, 0x54, 0x55, 0x56, 0x57, 0x58, 0x59, 0x5A, 0x63, 0x64, 0x65, 0x66, 0x67, 0x68, 0x69,
    0x6A, 0x73, 0x74, 0x75, 0x76, 0x77, 0x78, 0x79, 0x7A, 0x83, 0x84, 0x85, 0x86, 0x87, 0x88, 0x89,
    0x8A, 0x92, 0x93, 0x94, 0x95, 0x96, 0x97, 0x98, 0x99, 0x9A, 0xA2, 0xA3, 0xA4, 0xA5, 0xA6, 0xA7,
    0xA8, 0xA9, 0xAA, 0xB2, 0xB3, 0xB4, 0xB5, 0xB6, 0xB7, 0xB8, 0xB9, 0xBA, 0xC2, 0xC3, 0xC4, 0xC5,
    0xC6, 0xC7, 0xC8, 0xC9, 0xCA, 0xD2, 0xD3, 0xD4, 0xD5, 0xD6, 0xD7, 0xD8, 0xD9, 0xDA, 0xE1, 0xE2,
    0xE3, 0xE4, 0xE5, 0xE6, 0xE7, 0xE8, 0xE9, 0xEA, 0xF1, 0xF2, 0xF3, 0xF4, 0xF5, 0xF6, 0xF7, 0xF8,
    0xF9, 0xFA
};


// 标准 DC 色度霍夫曼表 (Bits, Values)
const uint8_t STD_DC_CHROM_BITS[17] = { 0, 0, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0 };
const uint8_t STD_DC_CHROM_VAL[12] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };

// 标准 AC 色度霍夫曼表 (Bits, Values)
const uint8_t STD_AC_CHROM_BITS[17] = { 0, 0, 2, 1, 2, 4, 4, 3, 4, 7, 5, 4, 4, 0, 1, 2, 119 };
const uint8_t STD_AC_CHROM_VAL[162] = { /* ... (省略具体值，与英文版相同) ... */
    0x00, 0x01, 0x02, 0x03, 0x11, 0x04, 0x05, 0x21, 0x31, 0x06, 0x12, 0x41, 0x51, 0x07, 0x61, 0x71,
    0x13, 0x22, 0x32, 0x81, 0x08, 0x14, 0x42, 0x91, 0xA1, 0xB1, 0xC1, 0x09, 0x23, 0x33, 0x52, 0xF0,
    0x15, 0x62, 0x72, 0xD1, 0x0A, 0x16, 0x24, 0x34, 0xE1, 0x25, 0xF1, 0x17, 0x18, 0x19, 0x1A, 0x26,
    0x27, 0x28, 0x29, 0x2A, 0x35, 0x36, 0x37, 0x38, 0x39, 0x3A, 0x43, 0x44, 0x45, 0x46, 0x47, 0x48,
    0x49, 0x4A, 0x53, 0x54, 0x55, 0x56, 0x57, 0x58, 0x59, 0x5A, 0x63, 0x64, 0x65, 0x66, 0x67, 0x68,
    0x69, 0x6A, 0x73, 0x74, 0x75, 0x76, 0x77, 0x78, 0x79, 0x7A, 0x82, 0x83, 0x84, 0x85, 0x86, 0x87,
    0x88, 0x89, 0x8A, 0x92, 0x93, 0x94, 0x95, 0x96, 0x97, 0x98, 0x99, 0x9A, 0xA2, 0xA3, 0xA4, 0xA5,
    0xA6, 0xA7, 0xA8, 0xA9, 0xAA, 0xB2, 0xB3, 0xB4, 0xB5, 0xB6, 0xB7, 0xB8, 0xB9, 0xBA, 0xC2, 0xC3,
    0xC4, 0xC5, 0xC6, 0xC7, 0xC8, 0xC9, 0xCA, 0xD2, 0xD3, 0xD4, 0xD5, 0xD6, 0xD7, 0xD8, 0xD9, 0xDA,
    0xE2, 0xE3, 0xE4, 0xE5, 0xE6, 0xE7, 0xE8, 0xE9, 0xEA, 0xF2, 0xF3, 0xF4, 0xF5, 0xF6, 0xF7, 0xF8,
    0xF9, 0xFA
};


// --- Z 字形扫描顺序 --- (与解码器相同)
const int G_ZIGZAG[BLOCK_SIZE * BLOCK_SIZE] = { /* ... (省略具体值，与英文版相同) ... */
     0,  1,  8, 16,  9,  2,  3, 10, 17, 24, 32, 25, 18, 11,  4,  5, 12, 19, 26, 33, 40, 48, 41, 34,
    27, 20, 13,  6,  7, 14, 21, 28, 35, 42, 49, 56, 57, 50, 43, 36, 29, 22, 15, 23, 30, 37, 44, 51,
    58, 59, 52, 45, 38, 31, 39, 46, 53, 60, 61, 54, 47, 55, 62, 63
};


// 预计算的 FDCT (正向离散余弦变换) 余弦值 (与 IDCT 数学上相同)
double G_FDCT_COS[BLOCK_SIZE][BLOCK_SIZE];

// --- 数据结构定义 ---

// 存储霍夫曼编码信息 (根据标准表预计算)
typedef struct {
    uint16_t code[256]; // 每个符号对应的霍夫曼码的数值
    uint8_t length[256]; // 每个符号对应的霍夫曼码的长度 (比特数)
} HuffmanCodeTable;

// 编码器使用的颜色分量信息
typedef struct {
    uint8_t id;             // 分量 ID (1=Y, 2=Cb, 3=Cr)
    uint8_t h_samp_factor; // 水平采样因子
    uint8_t v_samp_factor; // 垂直采样因子
    uint8_t quant_table_id; // 量化表 ID (0=亮度, 1=色度)
    int dc_pred;          // DC 系数预测值 (用于差分编码)
    // 霍夫曼表 ID (0=亮度, 1=色度) - 在这个基础编码器中根据分量 ID 隐式确定
} EncoderComponentInfo;

// 存储比特流写入状态的结构
typedef struct {
    FILE *fp;           // 指向输出 JPEG 文件
    uint32_t bit_buffer; // 比特缓冲器 (用 32 位存储可能较长的码字+数值)
    int bits_in_buffer; // 缓冲器中当前有多少比特是有效的
} BitstreamState;

// 编码器总信息结构
typedef struct {
    uint16_t width;     // 图像宽度
    uint16_t height;    // 图像高度
    uint8_t *rgb_data;  // 输入的 RGB 像素数据
    uint8_t *y_data;    // 亮度 (Y) 分量数据
    uint8_t *cb_data;   // 蓝色色度 (Cb) 分量数据 (二次采样后)
    uint8_t *cr_data;   // 红色色度 (Cr) 分量数据 (二次采样后)
    uint16_t cbcr_width; // 二次采样后色度平面的宽度
    uint16_t cbcr_height;// 二次采样后色度平面的高度

    EncoderComponentInfo components[MAX_COMPONENTS]; // 各分量信息
    int quality_factor; // 质量因子 (1-100)
    uint16_t quant_tables[2][BLOCK_SIZE * BLOCK_SIZE]; // 存储缩放后的量化表 (0: 亮度, 1: 色度)

    // 根据标准表预计算出的霍夫曼编码查找表
    HuffmanCodeTable dc_huff_codes[2]; // 0: 亮度, 1: 色度
    HuffmanCodeTable ac_huff_codes[2]; // 0: 亮度, 1: 色度

    BitstreamState bitstream; // 比特流写入状态

} JpegEncoderInfo;


// --- 工具函数 ---

// 将值限制在 0-255 范围内 (同解码器)
uint8_t clamp_byte(int val) {
    if (val < 0) return 0;
    if (val > 255) return 255;
    return (uint8_t)val;
}

// 向文件写入 1 个字节
void write_byte(uint8_t val, FILE *fp) {
    if (fputc(val, fp) == EOF) {
        perror("错误：向文件写入字节失败");
        exit(1);
    }
}

// 向文件写入 2 个字节 (大端序)
void write_word(uint16_t val, FILE *fp) {
    write_byte((uint8_t)((val >> 8) & 0xFF), fp); // 先写高位字节
    write_byte((uint8_t)(val & 0xFF), fp);      // 再写低位字节
}

// 写入一个 JPEG 标记 (0xFF 后跟标记代码)
void write_marker(uint8_t code, FILE *fp) {
    write_byte(0xFF, fp);
    write_byte(code, fp);
}

// 预计算 FDCT 余弦值 (与 IDCT 相同)
void init_fdct_cos() {
    for (int i = 0; i < BLOCK_SIZE; i++) {
        for (int j = 0; j < BLOCK_SIZE; j++) {
            G_FDCT_COS[i][j] = cos((2.0 * i + 1.0) * j * M_PI / (2.0 * BLOCK_SIZE));
        }
    }
}


// --- PPM 文件读取 ---
// 从 PPM P6 文件读取 RGB 图像数据
int read_ppm(const char *filename, JpegEncoderInfo *encoder) {
    FILE *fp_in = fopen(filename, "rb"); // 以二进制读取模式打开
    if (!fp_in) {
        perror("错误：打开输入 PPM 文件失败");
        return 0;
    }

    char magic[3]; // 存储 PPM 格式标识 ("P6")
    unsigned int width, height, maxval; // 存储宽度、高度、最大颜色值

    // 读取 PPM 文件头信息
    // %2s 读取 P6, %u 读取无符号整数, \n 跳过换行符
    if (fscanf(fp_in, "%2s\n%u %u\n%u\n", magic, &width, &height, &maxval) != 4) {
        fprintf(stderr, "错误：无效的 PPM 文件头格式。\n");
        fclose(fp_in);
        return 0;
    }

    // 检查是否是 P6 格式并且是 8 位图像
    if (magic[0] != 'P' || magic[1] != '6') {
        fprintf(stderr, "错误：输入文件不是 P6 格式的 PPM 文件。\n");
        fclose(fp_in);
        return 0;
    }
    if (maxval != 255) {
        fprintf(stderr, "错误：仅支持 8 位 PPM 文件 (最大颜色值为 255)。\n");
        fclose(fp_in);
        return 0;
    }

    // 存储图像尺寸
    encoder->width = (uint16_t)width;
    encoder->height = (uint16_t)height;
    size_t num_pixels = (size_t)width * height; // 总像素数
    size_t data_size = num_pixels * 3; // RGB 每个像素 3 字节

    printf("读取 PPM: %u x %u\n", width, height);

    // 分配内存存储 RGB 数据
    encoder->rgb_data = (uint8_t *)malloc(data_size);
    if (!encoder->rgb_data) {
        perror("错误：分配 RGB 数据内存失败");
        fclose(fp_in);
        return 0;
    }

    // 从文件中读取像素数据
    size_t read_count = fread(encoder->rgb_data, 1, data_size, fp_in);
    fclose(fp_in); // 读取完毕后关闭文件

    // 检查是否读取了完整的数据
    if (read_count != data_size) {
        fprintf(stderr, "错误：未能从 PPM 文件读取所有像素数据 (读取了 %zu / %zu 字节)。\n", read_count, data_size);
        free(encoder->rgb_data); // 释放已分配的内存
        encoder->rgb_data = NULL;
        return 0;
    }

    return 1; // 读取成功
}

// --- 颜色转换与二次采样 ---
// 将 RGB 数据转换为 YCbCr 并进行 4:2:0 色度二次采样
void rgb_to_ycbcr_and_subsample(JpegEncoderInfo *encoder) {
    size_t num_pixels = (size_t)encoder->width * encoder->height;

    // 分配 Y 分量内存 (与原图同尺寸)
    encoder->y_data = (uint8_t *)malloc(num_pixels);
    if (!encoder->y_data) { perror("分配 Y 分量内存失败"); exit(1); }

    // 计算 4:2:0 二次采样后 Cb/Cr 平面的尺寸 (宽高各减半，向上取整)
    encoder->cbcr_width = (encoder->width + 1) / 2;
    encoder->cbcr_height = (encoder->height + 1) / 2;
    size_t cbcr_num_pixels = (size_t)encoder->cbcr_width * encoder->cbcr_height;

    // 分配 Cb 和 Cr 分量内存
    encoder->cb_data = (uint8_t *)malloc(cbcr_num_pixels);
    encoder->cr_data = (uint8_t *)malloc(cbcr_num_pixels);
    if (!encoder->cb_data || !encoder->cr_data) { perror("分配 CbCr 分量内存失败"); exit(1); }

    printf("转换到 YCbCr 并进行二次采样 (4:2:0)...\n");
    printf("色度平面尺寸: %u x %u\n", encoder->cbcr_width, encoder->cbcr_height);

    // 遍历每个像素进行转换
    for (uint16_t y = 0; y < encoder->height; ++y) {
        for (uint16_t x = 0; x < encoder->width; ++x) {
            size_t rgb_idx = (y * encoder->width + x) * 3; // 当前像素 RGB 数据在数组中的起始索引
            size_t y_idx = y * encoder->width + x;        // 当前像素 Y 数据在数组中的索引
            uint8_t R = encoder->rgb_data[rgb_idx + 0];
            uint8_t G = encoder->rgb_data[rgb_idx + 1];
            uint8_t B = encoder->rgb_data[rgb_idx + 2];

            // 使用 JPEG 标准的 RGB 到 YCbCr 转换公式
            double Y_f =    0.299   * R + 0.587   * G + 0.114   * B;
            double Cb_f = - 0.168736* R - 0.331264* G + 0.5     * B + 128.0; // Cb/Cr 需要 +128 偏移
            double Cr_f =   0.5     * R - 0.418688* G - 0.081312* B + 128.0;

            // 存储 Y 分量 (钳位到 0-255)
            encoder->y_data[y_idx] = clamp_byte((int)round(Y_f));

            // 对 Cb 和 Cr 进行二次采样 (只在 2x2 块的左上角像素处计算并存储)
            if ((y % 2 == 0) && (x % 2 == 0)) {
                size_t cbcr_idx = (y / 2) * encoder->cbcr_width + (x / 2); // 计算在 Cb/Cr 平面中的索引

                // 对 2x2 区域的 Cb/Cr 值进行简单平均 (也可以只取左上角像素的值)
                // 这里采用平均法，考虑边界情况
                double cb_sum = 0, cr_sum = 0;
                int count = 0; // 计算实际参与平均的像素数
                for(int dy = 0; dy < 2; ++dy) {
                    for(int dx = 0; dx < 2; ++dx) {
                        uint16_t ny = y + dy; // 邻居像素的 y 坐标
                        uint16_t nx = x + dx; // 邻居像素的 x 坐标
                        // 确保邻居像素在图像范围内
                        if (ny < encoder->height && nx < encoder->width) {
                             size_t n_rgb_idx = (ny * encoder->width + nx) * 3;
                             uint8_t nR = encoder->rgb_data[n_rgb_idx + 0];
                             uint8_t nG = encoder->rgb_data[n_rgb_idx + 1];
                             uint8_t nB = encoder->rgb_data[n_rgb_idx + 2];
                             cb_sum += -0.168736* nR - 0.331264* nG + 0.5     * nB + 128.0;
                             cr_sum +=  0.5     * nR - 0.418688* nG - 0.081312* nB + 128.0;
                             count++;
                        }
                    }
                }
                // 存储平均后的 Cb/Cr 值 (钳位到 0-255)
                encoder->cb_data[cbcr_idx] = clamp_byte((int)round(cb_sum / count));
                encoder->cr_data[cbcr_idx] = clamp_byte((int)round(cr_sum / count));
            }
        }
    }
    // YCbCr 数据已生成，释放原始 RGB 数据内存
    free(encoder->rgb_data);
    encoder->rgb_data = NULL;
}


// --- 正向 DCT (FDCT) ---
// 对输入的 8x8 像素块执行 FDCT，输出 8x8 的 DCT 系数
void perform_fdct(const uint8_t *input_block, int16_t *output_coeffs) {
    double block[BLOCK_SIZE][BLOCK_SIZE]; // 用于存储 level-shifted 像素值
    double temp[BLOCK_SIZE][BLOCK_SIZE];  // 中间计算结果

    // 1. Level Shift: 将 0-255 的像素值减去 128，转换到 -128 到 127 的范围
    for (int y = 0; y < BLOCK_SIZE; y++) {
        for (int x = 0; x < BLOCK_SIZE; x++) {
            block[y][x] = (double)input_block[y * BLOCK_SIZE + x] - 128.0;
        }
    }

    // 2. 执行两次一维 FDCT (先对行，再对列)
    // 2a. 对每一行执行一维 FDCT
    for (int y = 0; y < BLOCK_SIZE; y++) {
        for (int u = 0; u < BLOCK_SIZE; u++) { // u 是频率索引
            double sum = 0.0;
            for (int x = 0; x < BLOCK_SIZE; x++) { // x 是空间索引
                sum += block[y][x] * G_FDCT_COS[x][u]; // 使用预计算的余弦值
            }
            double c_u = (u == 0) ? 1.0 / sqrt(2.0) : 1.0; // DCT-II 的归一化因子 C(u)
            temp[y][u] = c_u * sum / sqrt(BLOCK_SIZE / 2.0); // 标准 DCT-II 的归一化分母 sqrt(N/2)
        }
    }

    // 2b. 对临时结果的每一列执行一维 FDCT
    for (int u = 0; u < BLOCK_SIZE; u++) { // u 现在是列索引 (频率)
        for (int v = 0; v < BLOCK_SIZE; v++) { // v 是行索引 (频率)
            double sum = 0.0;
            for (int y = 0; y < BLOCK_SIZE; y++) { // y 是行索引 (空间)
                sum += temp[y][u] * G_FDCT_COS[y][v]; // 使用预计算的余弦值
            }
            double c_v = (v == 0) ? 1.0 / sqrt(2.0) : 1.0; // DCT-II 的归一化因子 C(v)
            // 将最终的 DCT 系数四舍五入并存储为 16 位整数
            output_coeffs[v * BLOCK_SIZE + u] = (int16_t)round(c_v * sum / sqrt(BLOCK_SIZE / 2.0));
        }
    }
}

// --- 量化 ---
// 对 8x8 的 DCT 系数块进行量化
void quantize(int16_t *coeffs, const uint16_t *quant_table) {
    for (int i = 0; i < BLOCK_SIZE * BLOCK_SIZE; ++i) {
        // 量化：用 DCT 系数除以量化表对应位置的值，并进行四舍五入取整
        double div = (double)coeffs[i] / quant_table[i];
        coeffs[i] = (int16_t)floor(div + 0.5); // 简单的四舍五入方法
        // 注意：JPEG 标准对四舍五入有具体规定，但这里简化处理
    }
}

// --- 霍夫曼编码辅助函数 ---

// 计算一个系数值需要多少比特来表示 (即它的 Category/Size)
int get_category(int value) {
    if (value == 0) return 0; // 0 不需要额外比特
    value = abs(value); // 取绝对值
    int category = 0;
    // 找到表示这个绝对值所需的最小比特数
    while ((1 << category) <= value) {
        category++;
    }
    return category;
}

// 获取系数的变长编码 (Variable-Length Integer, VLI)
// 正数 V = V; 负数 V = (1 << category) + V - 1
uint16_t get_vl_code(int value, int category) {
    if (value > 0) {
        return (uint16_t)value; // 正数直接返回
    } else {
        // 负数的编码方式类似二进制补码，但稍微不同
        return (uint16_t)((1 << category) + value - 1);
    }
}

// 根据标准霍夫曼表的 Bits 和 Values 数组，构建编码查找表 (符号 -> 码字+长度)
void build_huffman_code_tables(const uint8_t *bits, const uint8_t *values, HuffmanCodeTable *table) {
    memset(table, 0, sizeof(HuffmanCodeTable)); // 清零查找表

    uint16_t current_code = 0; // 当前生成的霍夫曼码的数值，从 0 开始
    int huff_val_index = 0;   // values 数组的当前索引

    // 遍历所有可能的码字长度 (1 到 16)
    for (int len = 1; len <= 16; ++len) {
        int count = bits[len]; // 当前长度 len 有多少个码字
        // 为当前长度的所有符号生成码字
        for (int i = 0; i < count; ++i) {
            uint8_t symbol = values[huff_val_index++]; // 获取当前码字对应的符号值
            if (symbol >= 256) { // 符号值不应超过 255
                 fprintf(stderr, "错误：霍夫曼符号值 >= 256 (%d)\n", symbol);
                 continue;
            }
            // 在查找表中记录该符号的码字数值和长度
            table->code[symbol] = current_code++;
            table->length[symbol] = len;
        }
        // 准备下一个长度的起始码字：(当前长度最后一个码字 + 1) << 1
        current_code <<= 1;
    }
}

// --- 比特流写入 ---

// 将指定数量的比特写入输出缓冲/文件
void write_bits(BitstreamState *bs, uint32_t bits, int n_bits) {
    if (n_bits == 0) return; // 不写入 0 比特

    // 将要写入的 n_bits 添加到缓冲区的末尾 (低位)
    // 1. 创建掩码，只保留 bits 的低 n_bits 位
    uint32_t mask = (1UL << n_bits) - 1;
    // 2. 将掩码后的比特左移，放到缓冲区的空闲位置
    bs->bit_buffer |= (bits & mask) << (32 - bs->bits_in_buffer - n_bits);
    // 3. 更新缓冲区中有效比特的数量
    bs->bits_in_buffer += n_bits;

    // 当缓冲区中的比特数 >= 8 时，将高位的字节写入文件
    while (bs->bits_in_buffer >= 8) {
        uint8_t byte_to_write = (uint8_t)(bs->bit_buffer >> 24); // 取出最高位的字节
        write_byte(byte_to_write, bs->fp); // 写入文件

        // 处理 Marker Stuffing：如果写入的字节是 0xFF，则在后面补充一个 0x00
        if (byte_to_write == 0xFF) {
            write_byte(0x00, bs->fp);
        }

        // 将缓冲区内容左移 8 位，并更新有效比特数
        bs->bit_buffer <<= 8;
        bs->bits_in_buffer -= 8;
    }
}

// 在编码结束时，将缓冲区中剩余的不足一个字节的比特写入文件
void flush_bits(BitstreamState *bs) {
    if (bs->bits_in_buffer > 0) {
        // JPEG 规范要求用 1 来填充剩余的比特位
        uint32_t pad_mask = (1UL << (8 - bs->bits_in_buffer)) - 1; // 创建 (8 - n) 个 1
        write_bits(bs, pad_mask, 8 - bs->bits_in_buffer); // 写入填充比特，这会触发最后一个字节的写入
    }
    // 理论上，调用 write_bits 后，如果 bits_in_buffer >= 8 会写入字节，最后应该为 0
    // 添加一个检查和强制清零，确保万无一失
    if (bs->bits_in_buffer != 0) {
        // 如果还有比特（理论上不应发生在此处，write_bits应已处理），强制写入
         uint8_t last_byte = (uint8_t)(bs->bit_buffer >> 24);
         write_byte(last_byte, bs->fp);
         if (last_byte == 0xFF) write_byte(0x00, bs->fp); // 别忘了 stuffing
         bs->bits_in_buffer = 0;
         bs->bit_buffer = 0;
    }
}


// --- JPEG 标记写入 ---

// 写入 APP0 (JFIF) 标记段
void write_APP0(FILE *fp) {
    write_marker(0xE0, fp); // APP0 标记代码
    write_word(16, fp);     // 段长度 = 16 字节 (固定)
    // 写入 "JFIF\0" 标识符
    write_byte('J', fp); write_byte('F', fp); write_byte('I', fp); write_byte('F', fp); write_byte(0, fp);
    // 写入 JFIF 版本号 (主版本 1, 次版本 1)
    write_byte(1, fp);
    write_byte(1, fp);
    // 写入密度单位 (0 = 无单位，仅宽高比)
    write_byte(0, fp);
    // 写入 X, Y 密度 (设为 1)
    write_word(1, fp);
    write_word(1, fp);
    // 写入缩略图宽高 (设为 0)
    write_byte(0, fp);
    write_byte(0, fp);
}

// 写入 DQT (定义量化表) 标记段
void write_DQT(JpegEncoderInfo *encoder) {
    write_marker(0xDB, encoder->bitstream.fp); // DQT 标记代码
    // 计算段长度: 2(长度字段) + 每个表(1字节信息 + 64字节数据)
    // 这里我们写入两个表 (亮度和色度)
    write_word(2 + (1 + 64) * 2, encoder->bitstream.fp);

    // 写入亮度量化表 (ID 0)
    write_byte(0x00, encoder->bitstream.fp); // 信息字节: 精度 0 (8位), 表 ID 0
    // 按 Z 字形顺序写入 64 个量化值
    for (int i = 0; i < 64; ++i) {
        write_byte((uint8_t)encoder->quant_tables[0][G_ZIGZAG[i]], encoder->bitstream.fp);
    }

    // 写入色度量化表 (ID 1)
    write_byte(0x01, encoder->bitstream.fp); // 信息字节: 精度 0 (8位), 表 ID 1
    // 按 Z 字形顺序写入 64 个量化值
    for (int i = 0; i < 64; ++i) {
        write_byte((uint8_t)encoder->quant_tables[1][G_ZIGZAG[i]], encoder->bitstream.fp);
    }
}

// 写入 SOF0 (帧开始 - 基线 DCT) 标记段
void write_SOF0(JpegEncoderInfo *encoder) {
    write_marker(0xC0, encoder->bitstream.fp); // SOF0 标记代码
    // 计算段长度: 2(长度) + 1(精度) + 2(高) + 2(宽) + 1(分量数) + 每个分量 3 字节信息
    write_word(2 + 1 + 2 + 2 + 1 + 3 * 3, encoder->bitstream.fp); // 假设固定为 3 分量
    write_byte(8, encoder->bitstream.fp);       // 图像精度 = 8 位
    write_word(encoder->height, encoder->bitstream.fp); // 图像高度
    write_word(encoder->width, encoder->bitstream.fp);  // 图像宽度
    write_byte(3, encoder->bitstream.fp);       // 颜色分量数 = 3 (Y, Cb, Cr)

    // 写入 Y 分量信息 (ID=1)
    write_byte(1, encoder->bitstream.fp);       // 分量 ID = 1
    write_byte(0x22, encoder->bitstream.fp);    // 采样因子: 水平=2, 垂直=2 (对应 4:2:0 的 Y)
    write_byte(0, encoder->bitstream.fp);       // 使用的量化表 ID = 0 (亮度)

    // 写入 Cb 分量信息 (ID=2)
    write_byte(2, encoder->bitstream.fp);       // 分量 ID = 2
    write_byte(0x11, encoder->bitstream.fp);    // 采样因子: 水平=1, 垂直=1 (对应 4:2:0 的 Cb)
    write_byte(1, encoder->bitstream.fp);       // 使用的量化表 ID = 1 (色度)

    // 写入 Cr 分量信息 (ID=3)
    write_byte(3, encoder->bitstream.fp);       // 分量 ID = 3
    write_byte(0x11, encoder->bitstream.fp);    // 采样因子: 水平=1, 垂直=1 (对应 4:2:0 的 Cr)
    write_byte(1, encoder->bitstream.fp);       // 使用的量化表 ID = 1 (色度)
}

// 写入 DHT (定义霍夫曼表) 标记段
void write_DHT(JpegEncoderInfo *encoder) {
    // 先计算 DHT 段的总长度
    // 长度 = 2(长度字段) + 4 * (1字节信息 + 16字节bits数组 + N字节values数组)
    int dc_lum_nvals = 0; for(int i=1; i<=16; ++i) dc_lum_nvals += STD_DC_LUM_BITS[i];
    int ac_lum_nvals = 0; for(int i=1; i<=16; ++i) ac_lum_nvals += STD_AC_LUM_BITS[i];
    int dc_chr_nvals = 0; for(int i=1; i<=16; ++i) dc_chr_nvals += STD_DC_CHROM_BITS[i];
    int ac_chr_nvals = 0; for(int i=1; i<=16; ++i) ac_chr_nvals += STD_AC_CHROM_BITS[i];
    uint16_t length = 2 + (1 + 16 + dc_lum_nvals) + (1 + 16 + ac_lum_nvals) +
                      (1 + 16 + dc_chr_nvals) + (1 + 16 + ac_chr_nvals);

    write_marker(0xC4, encoder->bitstream.fp); // DHT 标记代码
    write_word(length, encoder->bitstream.fp); // 写入总长度

    // --- 写入 DC 亮度霍夫曼表 (类型 0, ID 0) ---
    write_byte(0x00, encoder->bitstream.fp); // 信息字节: 类型=0 (DC), 表 ID=0
    for (int i = 1; i <= 16; ++i) write_byte(STD_DC_LUM_BITS[i], encoder->bitstream.fp); // 写入 bits 数组
    for (int i = 0; i < dc_lum_nvals; ++i) write_byte(STD_DC_LUM_VAL[i], encoder->bitstream.fp); // 写入 values 数组

    // --- 写入 AC 亮度霍夫曼表 (类型 1, ID 0) ---
    write_byte(0x10, encoder->bitstream.fp); // 信息字节: 类型=1 (AC), 表 ID=0
    for (int i = 1; i <= 16; ++i) write_byte(STD_AC_LUM_BITS[i], encoder->bitstream.fp);
    for (int i = 0; i < ac_lum_nvals; ++i) write_byte(STD_AC_LUM_VAL[i], encoder->bitstream.fp);

    // --- 写入 DC 色度霍夫曼表 (类型 0, ID 1) ---
    write_byte(0x01, encoder->bitstream.fp); // 信息字节: 类型=0 (DC), 表 ID=1
    for (int i = 1; i <= 16; ++i) write_byte(STD_DC_CHROM_BITS[i], encoder->bitstream.fp);
    for (int i = 0; i < dc_chr_nvals; ++i) write_byte(STD_DC_CHROM_VAL[i], encoder->bitstream.fp);

    // --- 写入 AC 色度霍夫曼表 (类型 1, ID 1) ---
    write_byte(0x11, encoder->bitstream.fp); // 信息字节: 类型=1 (AC), 表 ID=1
    for (int i = 1; i <= 16; ++i) write_byte(STD_AC_CHROM_BITS[i], encoder->bitstream.fp);
    for (int i = 0; i < ac_chr_nvals; ++i) write_byte(STD_AC_CHROM_VAL[i], encoder->bitstream.fp);
}

// 写入 SOS (扫描开始) 标记段
void write_SOS(JpegEncoderInfo *encoder) {
    write_marker(0xDA, encoder->bitstream.fp); // SOS 标记代码
    // 计算段长度: 2(长度) + 1(扫描分量数) + 每个分量 2 字节(ID+表ID) + 3 字节(谱选择)
    write_word(2 + 1 + 3 * 2 + 3, encoder->bitstream.fp);
    write_byte(3, encoder->bitstream.fp); // 本次扫描包含的分量数 = 3

    // 指定 Y 分量 (ID=1) 使用的霍夫曼表
    write_byte(1, encoder->bitstream.fp); // 分量选择器 = 1 (Y)
    write_byte(0x00, encoder->bitstream.fp); // DC 表 ID=0, AC 表 ID=0

    // 指定 Cb 分量 (ID=2) 使用的霍夫曼表
    write_byte(2, encoder->bitstream.fp); // 分量选择器 = 2 (Cb)
    write_byte(0x11, encoder->bitstream.fp); // DC 表 ID=1, AC 表 ID=1

    // 指定 Cr 分量 (ID=3) 使用的霍夫曼表
    write_byte(3, encoder->bitstream.fp); // 分量选择器 = 3 (Cr)
    write_byte(0x11, encoder->bitstream.fp); // DC 表 ID=1, AC 表 ID=1

    // 写入谱选择参数 (对于基线 DCT，通常是固定的)
    write_byte(0, encoder->bitstream.fp);   // 谱选择开始 (Ss) = 0 (DC)
    write_byte(63, encoder->bitstream.fp);  // 谱选择结束 (Se) = 63 (最后一个 AC 系数)
    write_byte(0, encoder->bitstream.fp);   // 逐次逼近位位置高位(Ah)=0, 低位(Al)=0
}

// --- 主编码循环 ---
// 对整个图像进行编码 (FDCT, 量化, 霍夫曼编码)
void encode_image(JpegEncoderInfo *encoder) {
    uint8_t input_block[BLOCK_SIZE * BLOCK_SIZE];   // 存储从 Y/Cb/Cr 平面提取的 8x8 像素块
    int16_t dct_coeffs[BLOCK_SIZE * BLOCK_SIZE];    // 存储 FDCT 后的系数
    int16_t quant_coeffs[BLOCK_SIZE * BLOCK_SIZE];  // 存储量化后的系数

    // 初始化各分量的 DC 系数预测值
    encoder->components[0].dc_pred = 0; // Y
    encoder->components[1].dc_pred = 0; // Cb
    encoder->components[2].dc_pred = 0; // Cr

    // 初始化比特流写入状态
    encoder->bitstream.bit_buffer = 0;
    encoder->bitstream.bits_in_buffer = 0;

    // 计算 MCU (最小编码单元) 的尺寸和数量 (基于 4:2:0 采样)
    // Y 分量的采样因子是 2x2，决定了 MCU 的像素尺寸
    uint16_t mcu_w_pixels = BLOCK_SIZE * encoder->components[0].h_samp_factor; // = 8 * 2 = 16
    uint16_t mcu_h_pixels = BLOCK_SIZE * encoder->components[0].v_samp_factor; // = 8 * 2 = 16
    // 计算图像按 MCU 划分的行列数
    uint16_t mcus_per_row = (encoder->width + mcu_w_pixels - 1) / mcu_w_pixels;
    uint16_t mcus_per_col = (encoder->height + mcu_h_pixels - 1) / mcu_h_pixels;

    printf("开始编码 MCUs (%u x %u)...\n", mcus_per_row, mcus_per_col);

    // 按 MCU 遍历图像
    for (uint16_t mcu_y = 0; mcu_y < mcus_per_col; ++mcu_y) {
        for (uint16_t mcu_x = 0; mcu_x < mcus_per_row; ++mcu_x) {

            // --- 处理 Y 分量的 4 个 8x8 块 (因为 Y 是 2x2 采样) ---
            EncoderComponentInfo *y_comp = &encoder->components[0];
            HuffmanCodeTable *dc_table_y = &encoder->dc_huff_codes[0]; // 亮度 DC 霍夫曼表
            HuffmanCodeTable *ac_table_y = &encoder->ac_huff_codes[0]; // 亮度 AC 霍夫曼表
            const uint16_t *q_table_y = encoder->quant_tables[0];      // 亮度量化表

            for (int v_samp = 0; v_samp < y_comp->v_samp_factor; ++v_samp) { // Y 在垂直方向有 2 个块
                for (int h_samp = 0; h_samp < y_comp->h_samp_factor; ++h_samp) { // Y 在水平方向有 2 个块
                    // 1. 从 Y 数据平面提取当前 8x8 块
                    uint16_t block_y_start = mcu_y * mcu_h_pixels + v_samp * BLOCK_SIZE;
                    uint16_t block_x_start = mcu_x * mcu_w_pixels + h_samp * BLOCK_SIZE;
                    for(int y=0; y<BLOCK_SIZE; ++y) {
                        for(int x=0; x<BLOCK_SIZE; ++x) {
                             uint16_t src_y = block_y_start + y;
                             uint16_t src_x = block_x_start + x;
                             // 处理边界情况：如果超出图像范围，用中间值 (128) 填充
                             if (src_y < encoder->height && src_x < encoder->width) {
                                 input_block[y*BLOCK_SIZE + x] = encoder->y_data[src_y * encoder->width + src_x];
                             } else {
                                 input_block[y*BLOCK_SIZE + x] = 128; // 填充
                             }
                        }
                    }

                    static int print_cnt = 2;
                    if(print_cnt) {
                        uint8_t *in = input_block;
                        printf("\n");
                        printf("input 8x8 block:\n");
                        for(int i = 0; i < 8; i++) {
                            for(int j = 0; j < 8; j++) {
                                printf("%4d ", in[i*8+j]);
                            }
                            printf("\n");
                        }
                    }

                    // 2. 执行 FDCT
                    perform_fdct(input_block, dct_coeffs);

                    if(print_cnt) {
                        int16_t *in = dct_coeffs;
                        printf("\n");
                        printf("after dct:\n");
                        for(int i = 0; i < 8; i++) {
                            for(int j = 0; j < 8; j++) {
                                printf("%4d ", in[i*8+j]);
                            }
                            printf("\n");
                        }
                    }

                    // 3. 量化
                    memcpy(quant_coeffs, dct_coeffs, sizeof(quant_coeffs)); // 复制一份用于量化
                    quantize(quant_coeffs, q_table_y);

                    if(print_cnt) {
                        int16_t *in = quant_coeffs;
                        printf("\n");
                        printf("after quantize:\n");
                        for(int i = 0; i < 8; i++) {
                            for(int j = 0; j < 8; j++) {
                                printf("%4d ", in[i*8+j]);
                            }
                            printf("\n");
                        }
                    }
                    print_cnt--;
                    if(print_cnt < 0)print_cnt = 0;

                    // 4. 对量化后的系数进行霍夫曼编码
                    // 4a. 编码 DC 系数 (差分编码)
                    int dc_diff = quant_coeffs[0] - y_comp->dc_pred; // 计算与前一块 DC 值的差值
                    y_comp->dc_pred = quant_coeffs[0]; // 更新 DC 预测值
                    int dc_category = get_category(dc_diff); // 获取差值的 Category (Size)
                    // 写入 DC 系数 Category 的霍夫曼码
                    write_bits(&encoder->bitstream, dc_table_y->code[dc_category], dc_table_y->length[dc_category]);
                    // 如果 Category > 0，则写入差值的 VLI 编码
                    if (dc_category > 0) {
                        write_bits(&encoder->bitstream, get_vl_code(dc_diff, dc_category), dc_category);
                    }

                    // 4b. 编码 AC 系数 (行程编码 RLE + 霍夫曼)
                    int zero_run = 0; // 记录连续 0 的个数
                    for (int k = 1; k < 64; ++k) { // 从第 1 个 AC 系数开始 (按 Z 字形顺序)
                        int coeff = quant_coeffs[G_ZIGZAG[k]]; // 获取 Z 字形扫描的下一个系数
                        if (coeff == 0) {
                            zero_run++; // 如果是 0，增加连续 0 的计数
                        } else {
                            // 遇到非零系数
                            // 先处理前面累计的 0：如果超过 15 个，需要写入 ZRL (16 个 0) 符号
                            while (zero_run >= 16) {
                                write_bits(&encoder->bitstream, ac_table_y->code[0xF0], ac_table_y->length[0xF0]); // 写入 ZRL (0xF0) 的霍夫曼码
                                zero_run -= 16;
                            }
                            // 获取非零系数的 Category (Size)
                            int ac_category = get_category(coeff);
                            // 构建 RLE 符号: 高 4 位是 run (连续 0 的个数)，低 4 位是 size (Category)
                            uint8_t rle_symbol = (zero_run << 4) | ac_category;
                            // 写入 RLE 符号的霍夫曼码
                            write_bits(&encoder->bitstream, ac_table_y->code[rle_symbol], ac_table_y->length[rle_symbol]);
                            // 写入非零系数的 VLI 编码
                            write_bits(&encoder->bitstream, get_vl_code(coeff, ac_category), ac_category);
                            zero_run = 0; // 重置连续 0 的计数
                        }
                    }
                    // 如果块的末尾有连续的 0，需要写入 EOB (End of Block) 符号
                    if (zero_run > 0) {
                        write_bits(&encoder->bitstream, ac_table_y->code[0x00], ac_table_y->length[0x00]); // 写入 EOB (0x00) 的霍夫曼码
                    }

                } // Y 水平块循环结束
            } // Y 垂直块循环结束


            // --- 处理 Cb 分量的 1 个 8x8 块 (因为 Cb 是 1x1 采样) ---
            EncoderComponentInfo *cb_comp = &encoder->components[1];
            HuffmanCodeTable *dc_table_c = &encoder->dc_huff_codes[1]; // 色度 DC 霍夫曼表
            HuffmanCodeTable *ac_table_c = &encoder->ac_huff_codes[1]; // 色度 AC 霍夫曼表
            const uint16_t *q_table_c = encoder->quant_tables[1];      // 色度量化表

            // 1. 从 Cb 数据平面提取 8x8 块 (注意使用 cbcr_width/height)
             uint16_t cb_block_y_start = mcu_y * BLOCK_SIZE; // MCU y 索引直接对应 Cb 块 y 索引
             uint16_t cb_block_x_start = mcu_x * BLOCK_SIZE; // MCU x 索引直接对应 Cb 块 x 索引
             for(int y=0; y<BLOCK_SIZE; ++y) {
                 for(int x=0; x<BLOCK_SIZE; ++x) {
                      uint16_t src_y = cb_block_y_start + y;
                      uint16_t src_x = cb_block_x_start + x;
                      if (src_y < encoder->cbcr_height && src_x < encoder->cbcr_width) {
                          input_block[y*BLOCK_SIZE + x] = encoder->cb_data[src_y * encoder->cbcr_width + src_x];
                      } else {
                          input_block[y*BLOCK_SIZE + x] = 128; // 边界填充
                      }
                 }
             }
             // 2. FDCT
             perform_fdct(input_block, dct_coeffs);
             // 3. 量化
             memcpy(quant_coeffs, dct_coeffs, sizeof(quant_coeffs));
             quantize(quant_coeffs, q_table_c);
             // 4. 霍夫曼编码 (DC)
             int dc_diff_cb = quant_coeffs[0] - cb_comp->dc_pred;
             cb_comp->dc_pred = quant_coeffs[0];
             int dc_category_cb = get_category(dc_diff_cb);
             write_bits(&encoder->bitstream, dc_table_c->code[dc_category_cb], dc_table_c->length[dc_category_cb]);
             if (dc_category_cb > 0) {
                 write_bits(&encoder->bitstream, get_vl_code(dc_diff_cb, dc_category_cb), dc_category_cb);
             }
             // 4b. 霍夫曼编码 (AC)
             int zero_run_cb = 0;
             for (int k = 1; k < 64; ++k) {
                 int coeff = quant_coeffs[G_ZIGZAG[k]];
                 if (coeff == 0) zero_run_cb++;
                 else {
                     while (zero_run_cb >= 16) { write_bits(&encoder->bitstream, ac_table_c->code[0xF0], ac_table_c->length[0xF0]); zero_run_cb -= 16; }
                     int ac_category = get_category(coeff);
                     uint8_t rle_symbol = (zero_run_cb << 4) | ac_category;
                     write_bits(&encoder->bitstream, ac_table_c->code[rle_symbol], ac_table_c->length[rle_symbol]);
                     write_bits(&encoder->bitstream, get_vl_code(coeff, ac_category), ac_category);
                     zero_run_cb = 0;
                 }
             }
             if (zero_run_cb > 0) write_bits(&encoder->bitstream, ac_table_c->code[0x00], ac_table_c->length[0x00]);


            // --- 处理 Cr 分量的 1 个 8x8 块 (因为 Cr 是 1x1 采样) ---
            EncoderComponentInfo *cr_comp = &encoder->components[2];
            // 在这个基础编码器中，Cr 使用与 Cb 相同的霍夫曼表和量化表
            // 1. 从 Cr 数据平面提取 8x8 块
             uint16_t cr_block_y_start = mcu_y * BLOCK_SIZE;
             uint16_t cr_block_x_start = mcu_x * BLOCK_SIZE;
              for(int y=0; y<BLOCK_SIZE; ++y) {
                 for(int x=0; x<BLOCK_SIZE; ++x) {
                      uint16_t src_y = cr_block_y_start + y;
                      uint16_t src_x = cr_block_x_start + x;
                      if (src_y < encoder->cbcr_height && src_x < encoder->cbcr_width) {
                          input_block[y*BLOCK_SIZE + x] = encoder->cr_data[src_y * encoder->cbcr_width + src_x];
                      } else {
                          input_block[y*BLOCK_SIZE + x] = 128; // 边界填充
                      }
                 }
             }
             // 2. FDCT
             perform_fdct(input_block, dct_coeffs);
             // 3. 量化 (使用色度量化表 q_table_c)
             memcpy(quant_coeffs, dct_coeffs, sizeof(quant_coeffs));
             quantize(quant_coeffs, q_table_c);
             // 4. 霍夫曼编码 (DC)
             int dc_diff_cr = quant_coeffs[0] - cr_comp->dc_pred;
             cr_comp->dc_pred = quant_coeffs[0];
             int dc_category_cr = get_category(dc_diff_cr);
             write_bits(&encoder->bitstream, dc_table_c->code[dc_category_cr], dc_table_c->length[dc_category_cr]);
             if (dc_category_cr > 0) {
                 write_bits(&encoder->bitstream, get_vl_code(dc_diff_cr, dc_category_cr), dc_category_cr);
             }
             // 4b. 霍夫曼编码 (AC)
             int zero_run_cr = 0;
             for (int k = 1; k < 64; ++k) {
                 int coeff = quant_coeffs[G_ZIGZAG[k]];
                 if (coeff == 0) zero_run_cr++;
                 else {
                     while (zero_run_cr >= 16) { write_bits(&encoder->bitstream, ac_table_c->code[0xF0], ac_table_c->length[0xF0]); zero_run_cr -= 16; }
                     int ac_category = get_category(coeff);
                     uint8_t rle_symbol = (zero_run_cr << 4) | ac_category;
                     write_bits(&encoder->bitstream, ac_table_c->code[rle_symbol], ac_table_c->length[rle_symbol]);
                     write_bits(&encoder->bitstream, get_vl_code(coeff, ac_category), ac_category);
                     zero_run_cr = 0;
                 }
             }
             if (zero_run_cr > 0) write_bits(&encoder->bitstream, ac_table_c->code[0x00], ac_table_c->length[0x00]);


        } // MCU 水平循环结束
         if ((mcu_y + 1) % 10 == 0 || mcu_y == mcus_per_col - 1) { // 每 10 行或最后一行打印进度
             printf("已编码 MCU 行 %u / %u\n", mcu_y + 1, mcus_per_col);
         }
    } // MCU 垂直循环结束

    // 编码完成后，清空比特流缓冲区中剩余的比特
    flush_bits(&encoder->bitstream);
    printf("完成扫描数据编码。\n");
}

// --- 量化表缩放 ---
// 根据质量因子缩放标准量化表
void scale_quant_table(const uint8_t *std_table, uint16_t *output_table, int quality_factor) {
    // JPEG 标准中推荐的质量因子到缩放因子的近似转换公式
    int scale_factor;
    if (quality_factor < 50) {
        scale_factor = 5000 / quality_factor;
    } else {
        scale_factor = 200 - 2 * quality_factor;
    }

    // 缩放并调整量化表中的每个值
    for (int i = 0; i < 64; ++i) {
        int temp = ((int)std_table[i] * scale_factor + 50) / 100; // 缩放并加 50 再除以 100 来实现四舍五入
        // 确保量化值在 1 到 255 之间 (JPEG 规定量化值不能为 0)
        if (temp <= 0) temp = 1;
        if (temp > 255) temp = 255; // 虽然是 uint16_t，但标准表缩放后通常不会超过 255
        output_table[i] = (uint16_t)temp;
    }
}


// --- 清理工作 ---
// 释放编码器分配的内存
void cleanup_encoder(JpegEncoderInfo *encoder) {
    // rgb_data 应该在转换后就被释放了，但以防万一检查一下
    if (encoder->rgb_data) free(encoder->rgb_data);
    if (encoder->y_data) free(encoder->y_data);
    if (encoder->cb_data) free(encoder->cb_data);
    if (encoder->cr_data) free(encoder->cr_data);
    // 注意：如果 bitstream.fp 已打开，应在调用此函数前关闭
}

// --- 主函数 ---
int main(int argc, char *argv[]) {
    // 检查命令行参数是否正确
    if (argc != 4) {
        fprintf(stderr, "用法: %s <输入.ppm> <输出.jpg> <质量(1-100)>\n", argv[0]);
        return 1;
    }

    // 初始化编码器信息结构
    JpegEncoderInfo encoder = {0};
    // 获取质量因子参数
    encoder.quality_factor = atoi(argv[3]);
    if (encoder.quality_factor < 1 || encoder.quality_factor > 100) {
        fprintf(stderr, "错误：质量因子必须在 1 到 100 之间。\n");
        return 1;
    }

    // 1. 读取输入的 PPM 图像文件
    if (!read_ppm(argv[1], &encoder)) {
        return 1; // 读取失败则退出
    }

    // 2. 初始化编码器所需的各种表
    init_fdct_cos(); // 预计算 FDCT 余弦值

    printf("根据质量因子 %d 缩放量化表...\n", encoder.quality_factor);
    scale_quant_table(STD_LUM_QUANT_TBL, encoder.quant_tables[0], encoder.quality_factor); // 缩放亮度量化表
    scale_quant_table(STD_CHROM_QUANT_TBL, encoder.quant_tables[1], encoder.quality_factor); // 缩放色度量化表

    printf("构建标准霍夫曼编码查找表...\n");
    // 构建亮度 DC/AC 和色度 DC/AC 的霍夫曼编码查找表
    build_huffman_code_tables(STD_DC_LUM_BITS, STD_DC_LUM_VAL, &encoder.dc_huff_codes[0]);
    build_huffman_code_tables(STD_AC_LUM_BITS, STD_AC_LUM_VAL, &encoder.ac_huff_codes[0]);
    build_huffman_code_tables(STD_DC_CHROM_BITS, STD_DC_CHROM_VAL, &encoder.dc_huff_codes[1]);
    build_huffman_code_tables(STD_AC_CHROM_BITS, STD_AC_CHROM_VAL, &encoder.ac_huff_codes[1]);

    // 3. 进行颜色空间转换 (RGB -> YCbCr) 和色度二次采样 (4:2:0)
    rgb_to_ycbcr_and_subsample(&encoder);

    // 4. 打开输出 JPEG 文件，准备写入比特流
    encoder.bitstream.fp = fopen(argv[2], "wb"); // 以二进制写入模式打开
    if (!encoder.bitstream.fp) {
        perror("错误：打开输出 JPEG 文件失败");
        cleanup_encoder(&encoder); // 清理已分配的内存
        return 1;
    }
    // 初始化比特流写入状态
    encoder.bitstream.bit_buffer = 0;
    encoder.bitstream.bits_in_buffer = 0;

    // 5. 写入 JPEG 文件头标记段
    write_marker(0xD8, encoder.bitstream.fp); // SOI (图像开始)
    write_APP0(encoder.bitstream.fp);         // APP0 (JFIF 头)
    write_DQT(&encoder);                      // DQT (量化表)
    write_SOF0(&encoder);                     // SOF0 (帧开始 - 基线)
    write_DHT(&encoder);                      // DHT (霍夫曼表)
    write_SOS(&encoder);                      // SOS (扫描开始)
                                              // SOS 之后紧跟着就是压缩的图像数据了

    // 6. 对图像数据进行编码 (核心处理循环)
    // 设置各分量信息 (ID, 采样因子, 量化表ID)
    encoder.components[0].id = 1; encoder.components[0].h_samp_factor = 2; encoder.components[0].v_samp_factor = 2; encoder.components[0].quant_table_id = 0; // Y
    encoder.components[1].id = 2; encoder.components[1].h_samp_factor = 1; encoder.components[1].v_samp_factor = 1; encoder.components[1].quant_table_id = 1; // Cb
    encoder.components[2].id = 3; encoder.components[2].h_samp_factor = 1; encoder.components[2].v_samp_factor = 1; encoder.components[2].quant_table_id = 1; // Cr
    encode_image(&encoder); // 执行编码

    // 7. 写入 EOI (图像结束) 标记
    write_marker(0xD9, encoder.bitstream.fp); // EOI

    // 8. 关闭输出文件并清理内存
    fclose(encoder.bitstream.fp); // 关闭文件
    cleanup_encoder(&encoder);   // 释放内存

    printf("成功将 %s 编码为 %s (质量: %d)\n", argv[1], argv[2], encoder.quality_factor);

    return 0; // 程序正常结束
}