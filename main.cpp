#include <iostream>
#include <string>
#include "ScoreFollowing/ScoreFollowing.h"

#include "wave.h"
#include "MFE/MFE.h"

using namespace std;
const int frame_len = 4096;
const int choose_template = 0;

int main()
{
    //读取.wav文件头信息
    // const char *path = "../data/new_example.wav";
    const char *path = "../data/new_example.wav";
    Wav wav;
    short *wavData = getWavMetaInfo(path, &wav, sizeof(wav));
    int wavLen = wav.data.Subchunk2Size / sizeof(short);
    int SampleRate = wav.fmt.SampleRate;
    int hopLen = 512;

    // 模板的最高频率是44100，如果使用不同的采样率必须重新训练得到模板
    // hopLen设置512和内部的乐谱跟随保持一致
    if (SampleRate != 44100 || hopLen != 512)
    {
        printf("the sampling rate must be 44100 and the hopLen must be 512");
        exit(1);
    }

    // 获取完整模板
    const char *templatePath = "../data/templateData/template11025.txt";
    double *templateData = getTemplateData(templatePath);

    // 挑选所需部分模板
    double *_templateData;
    // ScoreFollowing scoreFollowing;
    // vector<int> noteInScoreVector;
    // int *noteInScore;
    // if (choose_template)
    // {
    //     string scoreEventPath = "../data/scoreEvent/example.mat"; //由matlab readmidi_java.m得
    //     if (scoreFollowing.Init(scoreEventPath) == -1)
    //     {
    //         //读取乐谱信息文件失败
    //         cout << "读取乐谱信息文件失败" << endl;
    //         exit(1);
    //     }
    //     noteInScoreVector = scoreFollowing.GetNoteInScore();
    //     noteInScore = (int *)malloc(sizeof(int) * noteInScoreVector.size());
    //     for (int i = 0; i < noteInScoreVector.size(); ++i)
    //     {
    //         //            noteInScore[i] = noteInScoreVector[i];
    //         noteInScore[i] = noteInScoreVector[i] - 21; // noteInScoreVector为midi音高,对应MIDI值为21~108,转为0~87
    //     }
    //     int _newTemplateColumn;
    //     _templateData = chooseTemplate(templateData, noteInScore, noteInScoreVector.size(), &_newTemplateColumn);
    // }
    // else
    {
        _templateData = templateData; //  不挑选模板
        noteSize = 88;
    }

    int index = -frame_len / 2;
    int result_size = 88 * (wavLen / hopLen + 1);
    int(*pianoRoll88)[88] = (int(*)[88])malloc(sizeof(int) * result_size);
    double(*xH88)[88] = (double(*)[88])malloc(sizeof(double) * result_size);
    memset(pianoRoll88, 0, sizeof(int) * result_size);
    memset(xH88, 0, sizeof(double) * result_size);

    // 得到H
    for (;; index += hopLen)
    {
        // 一帧音频信息
        short samples[frame_len] = {0}; // 一帧
        // 前半帧 后半帧补零
        if (index < 0)
        {
            copy(wavData, &wavData[index + frame_len], &samples[-index]);
        }
        else if (index + frame_len < wavLen)
        {
            copy(&wavData[index], &wavData[index + frame_len], samples);
        }
        else if (index + frame_len / 2 < wavLen)
        {
            copy(&wavData[index], &wavData[wavLen], samples);
        }
        else
            break;

        double *xH = (double *)calloc(noteSize, sizeof(double));
        const char *c_windowPath = "../data/templateData/window4096_1024.txt";
        int *pianoRoll = MFEinterface(samples, wav.fmt.NumChannels, _templateData, xH, c_windowPath);

        // if (choose_template)
        // {
        //     for (int i = 0; i < noteSize; ++i)
        //     {
        //         pianoRoll88[(index + frame_len / 2) / hopLen][noteInScore[i] - 1] = pianoRoll[i];
        //         xH88[(index + frame_len / 2) / hopLen][noteInScore[i] - 1] = xH[i];
        //     }
        // }
        // else
        {
            for (int i = 0; i < 88; ++i)
            {
                pianoRoll88[(index + frame_len / 2) / hopLen][i] = pianoRoll[i];
                xH88[(index + frame_len / 2) / hopLen][i] = xH[i];
            }
        }
    }

    // cout << endl
    //      << endl
    //      << "pianoRoll | H :" << endl;
    // for (int i = 0; i < wavLen / hopLen + 1; i++)
    // {
    //     for (int j = 0; j < 88; ++j)
    //     {
    //         cout << pianoRoll88[i][j];
    //         cout << "|" << xH88[i][j] << "\t";
    //     }
    //     cout << i << endl;
    // }
    // 保存pianoRoll88和xH88
    FILE *pianoRoll88_file = fopen("../results/example_pianoRoll88.txt.", "w");
    FILE *xH88_file = fopen("../results/example_xH88.dat", "w");

    for (int j = 0; j < 88; j++)
    {
        for (int i = 0; i < wavLen / hopLen; i++)
        {
            // if (i != 0 && (i * hopLen * wav.fmt.NumChannels) % (SampleRate) < hopLen)
            //     fprintf(pianoRoll88_file, "|");
            // if (!pianoRoll88[i][j])
            // {
            //     fprintf(pianoRoll88_file, "_", pianoRoll88[i][j]);
            // }
            // else if (pianoRoll88[i][j] == 1)
            // {
            //     fprintf(pianoRoll88_file, "#", pianoRoll88[i][j]);
            // }
            // else
            // {
            //     fprintf(pianoRoll88_file, "%d ", pianoRoll88[i][j]);
            // }

            // if (xH88[i][j] >= 1e-6)
            // {
            //     fprintf(xH88_file, "%.4e\t", xH88[i][j]);
            // }
            // else
            // {
            //     fprintf(xH88_file, "%.4e\t", 0.0);
            // }

            fprintf(xH88_file, "%.4e ", xH88[i][j]);
        }
        // fprintf(pianoRoll88_file, "\n");
        fprintf(xH88_file, "\n");
    }
    free(wavData);
    // free(noteInScore);
    free(_templateData);
    free(pianoRoll88);
    free(xH88);
    return 0;
}
