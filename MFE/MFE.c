//
// Created by zhangqianyi on 2017/1/6.
//

#include "MFE.h"
#include "kiss_fft.h"
//#include <Accelerate/Accelerate.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define NULL ((void *)0)

extern double *resample(double *inputSignal, int upFactor, int downFactor, const char *path);


const int originNTemplate = 89; // 初始template模板的列数
int nTemplate = 89;        // template模板的列数
int noteSize = 88;           // 乐谱中音符个数，读取完乐谱之后确定
const int nPitch = 88;      // 音符个数
int maxPitches = 4;          // 一帧数据中最多出现的音符个数
const double beta = 0.6;    //nmf_beta算法中beta
const int niter = 15;    //nmf_beta算法中迭代次数
//double threshCoeff = 0.25;   // 阈值系数
double threshCoeff = 0.2;   // 阈值系数(真钢琴)

//double minThresh = 400;     // 阈值最小值
double minThresh = 240;     // 阈值最小值(真钢琴)

/**
 * @author zhangqianyi
 * @date   2016/12/6
 * @param  inputPath 表示模板文件路径
 * @return 返回保存模板数据的动态分配double型数组指针
 * @brief  读取模板数据
 */
// NSString *templatePath = [[NSBundle mainBundle] pathForResource:@"templatemaps" ofType:@"txt"];
// const char *path = [templatePath UTF8String];
// static double *templateData = getTemplateData(path);
double *getTemplateData(const char *inputPath) {
    int i;
    double *templateData = (double *) malloc(sizeof(double) * SPECTRUMLENGTH * nTemplate);   // 513*89 ? 2049*89
    FILE *templateFile;
    //从template.dat中读入频谱模板template (spectrumLength*nTemplate)
    if ((templateFile = fopen(inputPath, "rb")) == NULL) {
        printf("Cannot open template data file!\n");
        exit(1);
    }
    for (i = 0; i < SPECTRUMLENGTH * nTemplate; ++i) {
        fscanf(templateFile, "%lf", &templateData[i]);
    }
    fclose(templateFile);
    return templateData;
}

/**
 * @author zhangqianyi
 * @date 四月 2017
 *
 * @param templateData  完整的模板。要求：一套模板对应NPITCH个音符，列号为音符序号；多套模板拼接，各套对应的音符相同；若有静音模板，置于末尾，静音模板数 < NPITCH
 * @param noteInScore   乐谱包含的音符的序号
 * @param noteLength        noteInScore的长度
 * @param newTemplateColumn 返回模板的列数
 * @return newTemplate  乐谱包含的音符对应的模板。包含静音模板
 * @brief 根据乐谱包含的音符选择(抽取需要的)模板
 */
// static double *templateData = getTemplateData(path);
// vector<int> noteInScoreVector = scoreFollowing.GetNoteInScore();
// int *noteInScore = (int *)malloc(sizeof(int) * noteInScoreVector.size());
// for (int i = 0; i < noteInScoreVector.size(); ++i)	noteInScore[i] = noteInScoreVector[i];
// _templateData = chooseTemplate(templateData, noteInScore, noteInScoreVector.size(), &_newTemplateColumn);
double *chooseTemplate(double *templateData, int *noteInScore, int noteLength, int *newTemplateColumn) {
    int nTemplateSet = originNTemplate / nPitch;    // 有几套模板(初始template模板的列数89/音符个数88) 1
    int silenceTemplate = originNTemplate % nPitch; // 是否有静音模板(静音模板数) 1

    double *newTemplate;
    if (silenceTemplate == 0) { // 没有静音模板
        *newTemplateColumn = noteLength * nTemplateSet;
        newTemplate = (double *) malloc(sizeof(double) * SPECTRUMLENGTH * (*newTemplateColumn));
    } else { // 有静音模板
        *newTemplateColumn = noteLength * nTemplateSet + silenceTemplate;
        newTemplate = (double *) malloc(sizeof(double) * SPECTRUMLENGTH * (*newTemplateColumn));
    }

    int i, j, k;
    for (i = 0; i < nTemplateSet; ++i) {
        for (j = 0; j < noteLength; ++j) {
            int column = noteInScore[j] - 1;
            for (k = 0; k < SPECTRUMLENGTH; ++k) {
                // 第i个模板? 第j个音符 模板第k行
                newTemplate[i * noteLength + j + k * (*newTemplateColumn)]=templateData[i * nPitch + column + k * originNTemplate];
//                newTemplate[i * nTemplateSet + j + k * (*newTemplateColumn)] = templateData[i * nPitch + column + k * originNTemplate];
            }
        }
    }
    // 添加静音模板
    if (silenceTemplate != 0) {
        for (i = 0; i < SPECTRUMLENGTH; ++i) {
            newTemplate[*newTemplateColumn - 1 + i * (*newTemplateColumn)] = templateData[originNTemplate - 1 + i * originNTemplate];
        }
    }
    // 更新模板列数
    nTemplate = *newTemplateColumn;
    // 更新音符长度
    noteSize = noteLength;

    return newTemplate;
}

/**
* @author zhangqianyi
* @date 九月 2016
* @brief 多音调检测接口函数，输入原始录音数据，读入频谱模板处理之后调用多音调检测算法返回pianoRoll
* @param  short* originWav 输入原始录音数据，长度固定为4410，但是使用只取前4096个数据.
* @param  int wChannels 输入原始录音数据的声道数，默认为单声道数据.
* @param  path 采样窗路径
* @return int* pianoRoll 输出的结果.
*/
// short samples[4096] ≈ NSMutableArray *_MFEDataArray;//存放MFE需要的输入, 由多次录音数据拼接而成
// _templateData = chooseTemplate(templateData, noteInScore, noteInScoreVector.size(), &_newTemplateColumn);
// double *xH = (double *)calloc(noteSize, sizeof(double));
// NSString *windowPath = [[NSBundle mainBundle] pathForResource:@"window4096_1024" ofType:@"txt"];
// const char *c_windowPath = [windowPath UTF8String];
// int *pianoRoll = MFEinterface(samples, 1, _templateData, xH, c_windowPath);
int *MFEinterface(short *originWav, int wChannels, double *templateData, double *xH, const char *path) {
    int i;
    // 得到的原始录音数据，先归一化
    int framesize = 4096;
    double *normalizedData = (double *) malloc(framesize * sizeof(double));
    for (i = 0; i < framesize; ++i) {
        normalizedData[i] = (double) originWav[i] / 32767.0000; // 32767是16位整型数（int16）能表达的最大数.
    }

    // 读入频谱模板数据
    // double* templateData = (double *)malloc(sizeof(double)*SPECTRUMLENGTH*nTemplate);
    // getTemplateData(templateData);

    // 调用多音调检测算法

    int H0Flag = 0;
    double *H0 = (double *) calloc(nTemplate, sizeof(double));
    int *pianoRoll = (int *) calloc(noteSize, sizeof(int));
    // 对一帧数据normalizedData进行处理，得出一帧数据的pianoRoll值
    multiF0Estimation(normalizedData, templateData, H0Flag, H0, pianoRoll, xH, path);

    // 释放分配的内存
    free(normalizedData);
    normalizedData = NULL;
    free(H0);
    H0 = NULL;

    return pianoRoll;
}

/**
* @author zhangqianyi
* @date 七月 2016
* @brief 对一帧数据进行处理，得出一帧数据的pianoRoll值.
* @param  xFrame 输入的一帧数据.
* @param  templateData 输入的模板数据.
* @param  pianoRoll 输出的结果.
* @return 返回value.
*/
void multiF0Estimation(double *xFrame, double *templateData, int H0Flag, double *H0, int *pianoRoll, double *xH,
                       const char *path) {
    double *outframe = (double *) malloc(4096 * sizeof(double));
    outframe = resample(xFrame, 1024, 4096, path);
    int i, j;
// 先判断xFrame是否全为零，全为零就将pianoRoll置0返回
    double xFrameSum = 0.0;
    //for (i = 0; i < FRAMESIZE; ++i) xFrameSum += fabs(xFrame[i]);
    for (i = 0; i < FRAMESIZE; ++i) xFrameSum += fabs(outframe[i]);
    if (xFrameSum < 1e-8) {
        printf("xFrame is all zero!\n");
        for (i = 0; i < noteSize; ++i) {
            pianoRoll[i] = 0;
            xH[i] = 0;
        }
        return;
    }

    // 计算信号频谱
    double *spectrum = (double *) malloc(sizeof(double) * SPECTRUMLENGTH);
    spectrumOfFrame(outframe, spectrum);

    // 判断spectrum 中是否含零，包含零就将pianoRoll置0返回

    for (i = 0; i < SPECTRUMLENGTH; ++i) {
//        printf("Spectrum has zero!\n");
        if ((spectrum[i] - 0) < 1e-8) {
            for (j = 0; j < noteSize; ++j) {
                pianoRoll[j] = 0;
                xH[j] = 0;
            }
            // 先释放分配的内存，再返回，防止内存泄漏
            free(spectrum);
            spectrum = NULL;
            return;
        }
    }

    //频谱因式分解得到H
    double *H = (double *) calloc(nTemplate, sizeof(double)); // 定义H并初始化为0
    nmf_beta(spectrum, SPECTRUMLENGTH, nTemplate, templateData, H, H0Flag, H0, beta, niter);
    //nmf_beta_cblas(spectrum, SPECTRUMLENGTH, nTemplate, templateData, H, H0Flag, H0, beta, niter);


    //有多套模板时，将对应于同一音调的pitch activity相加
    for (i = 1; i < nTemplate / noteSize; ++i) {
        for (j = 0; j < noteSize; ++j) {
            H[j] += H[j + i * noteSize];
        }
    }

    for (i = 0; i < noteSize; ++i) xH[i] = H[i];
    // 得到H峰值maxH
    double maxH = H[0];
    for (i = 1; i < noteSize; ++i) {
        if (H[i] > maxH) {
            maxH = H[i];
        }
    }
    // 计算最小阈值
    double thisThresh = maxH * threshCoeff > minThresh ? maxH * threshCoeff : minThresh;

    // 加阈值   得到0 1形式的pianoRoll
    for (i = 0; i < noteSize; ++i) {
        pianoRoll[i] = H[i] > thisThresh;
    }

    // 限制一帧数据中出现的音符的个数，只取最大的前maxPitches个，其余的pianoRoll置0
/*
    int pianoRollSum = 0;
    for (i = 0; i < noteSize; ++i) {
        if (pianoRoll[i] == 1) ++pianoRollSum;
    }
    if (pianoRollSum > maxPitches) { // 先对H前notesize位排序(从大到小)，然后将下标maxPitches后面的pianoRoll设置为0
        for (i = 0; i < noteSize; ++i) pianoRoll[i] = 0;
        QSort(H, 0, noteSize - 1);
        int index = 0;
        for (index = 0; index < maxPitches; ++index) {
            for (i = 0; i < noteSize; ++i) {
                if (xH[i] == H[index]) {
                    pianoRoll[i] = 1;  // 将最大的前maxPitches位置为1
                    break;
                }
            }
        }
    }
*/
    free(H);H = NULL;
    free(spectrum);spectrum = NULL;
}

/**
 * @author zhangqianyi
 * @date 2016/7/12
 * @brief 对一帧信号进行FFT，计算振幅谱。 窗函数为Hamming窗，DFT偶对称；若帧长小于fft长度，信号加窗后末尾补0。
 * @param  xFrame 一帧信号（长度为帧长。若帧长为奇数，窗函数不满足DFT偶对称）
 * @param  spectrum 输出振幅谱地址（前半段, fftSize/2+1）
 */
void spectrumOfFrame(double *xFrame, double *spectrum) {
    int i;
    double *win = (double *) calloc(FRAMESIZE, sizeof(double));
    hammingWindow(win, FRAMESIZE);
    //加汉明窗
    double *xFft = (double *) calloc(FFTSIZE, sizeof(double));
    for (i = 0; i < FRAMESIZE; ++i) {
        xFft[i] = xFrame[i] * win[i];
    }
    if (FFTSIZE > FRAMESIZE) //加窗后超出帧的部分补0
    {
        for (i = FRAMESIZE; i < FFTSIZE; ++i) {
            xFft[i] = 0.0;
        }
    }
    free(win);
    win = NULL;

    //kiss-fft
    //input
    kiss_fft_cpx *kiss_xFft = (kiss_fft_cpx *) malloc(sizeof(kiss_fft_cpx) * FFTSIZE);
    for (i = 0; i < FFTSIZE; ++i) {
        kiss_xFft[i].r = xFft[i];
        kiss_xFft[i].i = 0.0;
    }
    free(xFft);
    xFft = NULL;

    //output
    kiss_fft_cpx *fftFrame = (kiss_fft_cpx *) malloc(sizeof(kiss_fft_cpx) * FFTSIZE);
    kiss_fft_cfg cfg = kiss_fft_alloc(FFTSIZE, 0, 0, 0); // 为fft分配所有必要的存储空间

    kiss_fft(cfg, kiss_xFft, fftFrame);

    //计算频谱振幅
    for (i = 0; i < SPECTRUMLENGTH; ++i) {
        spectrum[i] = sqrt(fftFrame[i].r * fftFrame[i].r + fftFrame[i].i * fftFrame[i].i); // 计算幅值
    }

    free(kiss_xFft);
    free(cfg);
    free(fftFrame);
    fftFrame = NULL;
}

/**
* @author zhangqianyi
* @date 七月 2016
* @brief 获得一个指定大小的hamming窗.
* @param  w 返回的窗函数值.
* @param  windowSize 窗函数的大小。
*/
void hammingWindow(double *w, int windowSize) {
    int i;
    double pi = 3.14159265358979323846;
    for (i = 0; i < windowSize; ++i) {
        w[i] = 0.54 - 0.46 * cos(2 * pi * i / (windowSize - 1));
    }
}

/**
* @author zhangqianyi
* @date 七月 2016
* @brief nmf_beta非负矩阵分解函数（data = W * H，已知data和W来得到H）.
* @param data：要分解的矩阵
* @param rows：W的行数
* @param columns：W的列数
* @param W：矩阵分解的频谱模板templateData
* @param H：矩阵分解的结果
* @param beta：用于矩阵分解的beta参数
* @param maxiter：循环计算的次数
* @param H0Flag：是否有H0输入，1为有H0，0为没有H0
* @param H0：H0数据
*/
// 	nmf_beta(spectrum, SPECTRUMLENGTH, nTemplate=89, templateData,     H,   H0Flag=0,         H0,        beta=0.6,niter=15);
void nmf_beta(double *data, int rows, int columns, double *W, double *H, int H0Flag, double *H0, double beta, int maxiter) {
    int i, j, k;
    // 初始化H
    if (H0Flag == 0) { // 没有H0输入，则随机初始化
        for (i = 0; i < columns; ++i) {
            // H[i] = 1;  // this is for test，the code below is for release
            H[i] = rand() / ((double) RAND_MAX + 1);
        }
    } else if (H0Flag == 1) { // 有H0输入，则将H初始化为H0
        for (i = 0; i < columns; ++i) {
            H[i] = H0[i];
        }
    }

        // 将模板标准化
//        double* mrowsum = (double*)calloc(columns, sizeof(double));
//        for (i = 0; i < columns; ++i) {
//            for (j = 0; j < rows; ++j) {
//                mrowsum[i] += W[j * columns + i]; // 每一列的总和
//            }
//        }
//        for (i = 0; i < rows; ++i) {
//            for (j = 0; j < columns; ++j) {
//                W[i * columns + j] /= mrowsum[j];
//            }
//        }
//        free(mrowsum); mrowsum = NULL;

    // 初始化R
    double *R = (double *) calloc(rows, sizeof(double));
    for (i = 0; i < rows; ++i) {
        for (k = 0; k < columns; ++k) {
            R[i] += W[i * columns + k] * H[k];  // 二维矩阵相乘
        }
    }

    double *temp1 = (double *) calloc(rows, sizeof(double));
    double *temp2 = (double *) calloc(rows, sizeof(double));
    double *temp3 = (double *) calloc(columns, sizeof(double));
    double *temp4 = (double *) calloc(columns, sizeof(double));

    // 迭代
    double myeps = 1.0e-20; // 精度
    int it; // 迭代次数
    for (it = 0; it < maxiter; ++it) {
        //update R
        memset(R, 0, sizeof(double) * rows); // 将R清零
        for (i = 0; i < rows; ++i) {
            for (k = 0; k < columns; ++k) {
                R[i] += W[i * columns + k] * H[k];
            }
        }
        //update H
        // 初始化中间变量，内存清零
        memset(temp1, 0, sizeof(double) * rows);
        memset(temp2, 0, sizeof(double) * rows);
        memset(temp3, 0, sizeof(double) * columns);
        memset(temp4, 0, sizeof(double) * columns);
        for (i = 0; i < rows; ++i) {
            temp1[i] = pow(R[i], beta - 2) * data[i];
            temp2[i] = pow(R[i], beta - 1);
        }
        for (k = 0; k < rows; ++k) {
            for (i = 0; i < columns; ++i) {
                temp3[i] += W[k * columns + i] * temp1[k];
                temp4[i] += W[k * columns + i] * temp2[k];
            }
        }
        for (i = 0; i < columns; ++i) {
            if (temp4[i] < myeps) {
                temp4[i] = myeps;
            }
            H[i] = H[i] * (temp3[i] / temp4[i]);
        }
    }

    free(temp1);temp1 = NULL;
    free(temp2);temp2 = NULL;
    free(temp3);temp3 = NULL;
    free(temp4);temp4 = NULL;
    free(R);R = NULL;
}


/**
* @author zhangqianyi
* @date   2017/4/24
* @param a     需要排序的数组
* @param low   初始化为0
* @param high  初始化为a.length-1
* @brief 快排，数组为double型
*/
void QSort(double *a, int low, int high) {
    if (low < high) {
        int i = low, j = high;
        double pivotkey = a[i];
        while (i < j) {
            while (i < j && a[j] < pivotkey) --j; // 找到大于等于p的a[j]
            a[i] = a[j];
            while (i < j && a[i] > pivotkey) ++i; // 找到小于等于p的a[i]
            a[j] = a[i];
        }
        a[i] = pivotkey;  // 将所有大于等于a[0]的数放在数组左边,小于等于a[0]的数放在右边,a[0]放在分界处

        QSort(a, low, i - 1);
        QSort(a, i + 1, high);
    }
}

#pragma mark - USE BLAS

//void nmf_beta_cblas(double *data, int rows, int columns, double *W, double *H, int H0Flag, double *H0, double beta,
//                    int maxiter) {
//    // clock_t a = clock();
//    int i, j, k;
//
//#pragma mark - H
//    // 初始化H
//    if (H0Flag == 0) { // 没有H0输入，则随机初始化
//        for (i = 0; i < columns; ++i) {
//            H[i] = 1;  // this is for test，the code below is for release
//            // H[i] = rand() / ((double)RAND_MAX + 1);
//        }
//    } else if (H0Flag == 1) { // 有H0输入，则将H初始化为H[0]
//        for (i = 0; i < columns; ++i) {
//            H[i] = H0[i];
//        }
//    }
//
//    // // 将W标准化 耗时3-4ms
//    // double* mRowSum = (double*)calloc(columns, sizeof(double));
//    // for (i = 0; i < columns; ++i) {
//    // 	for (j = 0; j < rows; ++j) {
//    // 		mRowSum[i] += W[j * columns + i]; // 每一列的总和
//    // 	}
//    // }
//    // for (i = 0; i < rows; ++i) {
//    // 	for (j = 0; j < columns; ++j) {
//    // 		W[i * columns + j] /= mRowSum[j];
//    // 	}
//    // }
//    // free(mRowSum); mRowSum = NULL;
//
//    // 初始化R
//    double *R = (double *) calloc(rows, sizeof(double));
//
//    //for (i = 0; i < rows; ++i) {
//    //	for (k = 0; k < columns; ++k) {
//    //		R[i] += W[i * columns + k] * H[k];  // 二维矩阵相乘
//    //	}
//    //}
//
//    // 使用cblas库，下面两行是介绍如何使用，每个参数的意义
//    // cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, A, k, B, n, beta, C, n);
//    // m=rows, n=1, k=columns, C=R, A=W, B=H, alpha=1, beta=0
//    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, rows, 1, columns, 1, W, columns, H, 1, 0, R, 1);
//
//
//    // printf("准备nmf耗时 = %d\t", (int)(clock() - a));
//    // a = clock();
//    double *temp1 = (double *) calloc(rows, sizeof(double));
//    double *temp2 = (double *) calloc(rows, sizeof(double));
//    double *temp3 = (double *) calloc(columns, sizeof(double));
//    double *temp4 = (double *) calloc(columns, sizeof(double));
//
//    // 迭代
//    double myeps = 1.0e-20; // 精度
//    int it; // 迭代次数
//    for (it = 0; it < maxiter; ++it) {
//        //update R
//        memset(R, 0, sizeof(double) * rows); // 将R清零
//
//        /*for (i = 0; i < rows; ++i) {
//         for (k = 0; k < columns; ++k) {
//         R[i] += W[i * columns + k] * H[k];
//         }
//         }*/
//        // 使用cblas库
//        // cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, A, k, B, n, beta, C, n);
//        // m=rows, n=1, k=columns, C=R, A=W, B=H, alpha=1, beta=0
//        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, rows, 1, columns, 1, W, columns, H, 1, 0, R, 1);
//
//        //update H
//        // 初始化中间变量，内存清零
//        memset(temp1, 0, sizeof(double) * rows);
//        memset(temp2, 0, sizeof(double) * rows);
//        memset(temp3, 0, sizeof(double) * columns);
//        memset(temp4, 0, sizeof(double) * columns);
//        for (i = 0; i < rows; ++i) {
//            temp1[i] = pow(R[i], beta - 2) * data[i];
//            temp2[i] = pow(R[i], beta - 1);
//        }
//
//        /*for (k = 0; k < rows; ++k) {
//         for (i = 0; i < columns; ++i) {
//         temp3[i] += W[k * columns + i] * temp1[k];
//         temp4[i] += W[k * columns + i] * temp2[k];
//         }
//         }*/
//        // cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, A, k, B, n, beta, C, n);
//        // m=1, n=columns, k=rows, C=temp3, A=temp1, B=W
//        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 1, columns, rows, 1, temp1, rows, W, columns, 0, temp3,
//                    columns);
//        // m=1, n=columns, k=rows, C=temp4, A=temp2, B=W
//        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 1, columns, rows, 1, temp2, rows, W, columns, 0, temp4,
//                    columns);
//
//        for (i = 0; i < columns; ++i) {
//            if (temp4[i] < myeps) {
//                temp4[i] = myeps;
//            }
//            H[i] = H[i] * (temp3[i] / temp4[i]);
//        }
//    }
//
//    free(temp1);
//    temp1 = NULL;
//    free(temp2);
//    temp2 = NULL;
//    free(temp3);
//    temp3 = NULL;
//    free(temp4);
//    temp4 = NULL;
//
//    // printf("nmf耗时 = %d\n", (int)(clock() - a));
//    free(R);
//    R = NULL;
//}

/**
 @brief BLAS计算矩阵相乘
 @param rowsA 输入矩阵的行
 @param columnsA 输入矩阵的列
 @param alpha 输入矩阵和输入向量的系数，一般为1
 @param matrixA 输入矩阵
 @param vector 输入向量
 @param beta 输出矩阵的偏移，一般为0
 @param matrexC 输出矩阵
 */
//void matrixMultiply(int rowsA, int columnsA, float alpha, float *matrixA, float *vector, float beta, float *matrexC) {
//    cblas_sgemv(CblasRowMajor, CblasNoTrans, rowsA, columnsA, alpha, matrixA, columnsA, vector, 1, beta, matrexC, 1);
//}
