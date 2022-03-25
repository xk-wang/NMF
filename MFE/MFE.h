#pragma once
#ifndef MFE_H_
#define MFE_H_
#ifdef __cplusplus
extern "C"
{
#endif

  extern const int originNTemplate; // 初始template模板的列数
  extern int nTemplate;             // template模板的列数
  extern int noteSize;              // 乐谱中音符个数，读取完乐谱之后确定
  extern int maxPitches;            // 一帧数据中最多出现的音符个数
  extern const int nPitch;          // 音符个数，这个是常数88

//#define FFTSIZE 4096                // fft长度
//#define SPECTRUMLENGTH 2049         // 谱长度（fftSize/2+1）
//#define HOPSIZE 441                 // 帧移
//#define FRAMESIZE 4096              // 帧长（samples）
#define FFTSIZE 1024       // fft长度
#define SPECTRUMLENGTH 513 // 谱长度（fftSize/2+1）
#define HOPSIZE 110        // 帧移
#define FRAMESIZE 1024     // 帧长（samples）

  extern const double beta; // nmf_beta算法中beta
  extern const int niter;   // nmf_beta算法中迭代次数

  //    extern double threshCoeff;    // 阈值系数
  //    extern double minThresh;      // 阈值最小值

  /**
   * @author zhangqianyi
   * @date 九月 2016
   * @brief 多音调检测接口函数，输入原始录音数据，读入频谱模板处理之后调用多音调检测算法返回pianoRoll
   * @param  short* originWav 输入原始录音数据，长度固定为4410，但是使用只取前4096个数据.
   * @param  int wChannels 输入原始录音数据的声道数，默认为单声道数据.
   * @return int* pianoRoll 输出的结果.
   */
  int *MFEinterface(short *originWav, int wChannels, double *templateData, double *xH, const char *path);

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
                         const char *path);

  /**
   * @author zhangqianyi
   * @date 2016/7/12
   * @brief 对一帧信号进行FFT，计算振幅谱。 窗函数为Hamming窗，DFT偶对称；若帧长小于fft长度，信号加窗后末尾补0。
   * @param  xFrame 一帧信号（长度为帧长。若帧长为奇数，窗函数不满足DFT偶对称）
   * @return 振幅谱地址（前半段, fftSize/2+1）
   */
  void spectrumOfFrame(double *xFrame, double *spectrum);

  /**
   * @author zhangqianyi
   * @date 七月 2016
   * @brief 获得一个指定大小的hamming窗.
   * @param  w 返回的窗函数值.
   * @param  windowSize 窗函数的大小。
   */
  void hammingWindow(double *w, int windowSize);

  /**
   * @author zhangqianyi
   * @date 七月 2016
   * @brief nmf_beta非负矩阵分解函数（data = W * H，已知data和W来得到H）.
   * @param data：要分解的矩阵
   * @param rows：W的行数
   * @param classes：W的列数
   * @param W：矩阵分解的频谱模板
   * @param H：矩阵分解的结果
   * @param beta：用于矩阵分解的beta参数
   * @param maxiter：循环计算的次数
   * @param H0Flag：是否有H0输入，1为有H0，0为没有H0
   * @param H0：H0数据
   */
  void
  nmf_beta(double *data, int rows, int column, double *W, double *H, int H0Flag, double *H0, double beta, int maxiter);

  /**
   * @author zhangqianyi
   * @date 四月 2017
   *
   * @param templateData  完整的模板。要求：一套模板对应NPITCH个音符，列号为音符序号；多套模板拼接，各套对应的音符相同；若有静音模板，置于末尾，静音模板数<NPITCH
   * @param noteInScore   乐谱包含的音符的序号
   * @param noteLength      noteInScore的长度
   * @param newTemplateColumn 返回模板的列数
   * @return newTemplate  乐谱包含的音符对应的模板。包含静音模板
   * @brief 根据乐谱包含的音符选择模板
   */
  double *chooseTemplate(double *templateData, int *noteInScore, int noteLength, int *newTemplateColumn);

  /**
   * @author zhangqianyi
   * @date   2016/12/6
   * @param  inputPath 表示模板文件路径
   * @return 返回保存模板数据的动态分配double型数组指针
   * @brief  读取模板数据
   */
  double *getTemplateData(const char *inputPath);

  void QSort(double *a, int low, int high);

  /**
  * @author zhangqianyi
  * @date 七月 2016
  * @brief nmf_beta非负矩阵分解函数（data = W * H，已知data和W来得到H）.
  * @param data：要分解的矩阵
  * @param rows：W的行数
  * @param classes：W的列数
  * @param W：矩阵分解的频谱模板
  * @param H：矩阵分解的结果
  * @param beta：用于矩阵分解的beta参数
  * @param maxiter：循环计算的次数


  尝试使用cblas库来完成矩阵相乘运算


  */

  void nmf_beta_cblas(double *data, int rows, int columnss, double *W, double *H, int H0Flag, double *H0, double beta,
                      int maxiter); //非负矩阵分解

  void matrixMultiply(int rowsA, int columnsA, float alpha, float *matrixA, float *vector, float beta, float *matrexC);
#ifdef __cplusplus
}
#endif

#endif // !MFE_H_
