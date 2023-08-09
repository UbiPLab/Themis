#include <iostream>
#include <opencv2/opencv.hpp>
#include <cstdlib>
#include <ctime>

using namespace std;
using namespace cv;

int main()
{
    int alpha = 2;
    int len = 2048;
    int wid = 2048;
    int m_num = 38;
    int b_size = m_num;
    int p_size = 2;
    int w_size = m_num * p_size;
    int macro_size = w_size * 2 * 2;
    int Lnum = len / macro_size;
    int Wnum = wid / macro_size;

    Mat I0 = imread("./input/example.jpg");
    cout << "---------loading---------\n";
    if (I0.empty() || !I0.data)
    {
        cout << "Loading input image failed!" << endl;
        return -1;
    }
    else
    {
        cout << "Loading input image succsessfully." << endl;
    }
    int tmp_height = I0.rows;
    int tmp_width = I0.cols;
    vector<vector<int>> tmp_r(tmp_height, vector<int>(tmp_width, 0));
    vector<vector<int>> tmp_g(tmp_height, vector<int>(tmp_width, 0));
    vector<vector<int>> tmp_b(tmp_height, vector<int>(tmp_width, 0));
    for (int i = 0; i < tmp_height; ++i)
    {
        for (int j = 0; j < tmp_width; ++j)
        {
            tmp_b[i][j] = I0.at<Vec3b>(i, j)[0];
            tmp_g[i][j] = I0.at<Vec3b>(i, j)[1];
            tmp_r[i][j] = I0.at<Vec3b>(i, j)[2];
        }
    }

    vector<vector<double>> tmp_Y(tmp_height, vector<double>(tmp_width, 0));
    vector<vector<double>> tmp_Cb(tmp_height, vector<double>(tmp_width, 0));
    vector<vector<double>> tmp_Cr(tmp_height, vector<double>(tmp_width, 0));
    for (int i = 0; i < tmp_height; ++i)
    {
        for (int j = 0; j < tmp_width; ++j)
        {
            tmp_Y[i][j] = 0.299 * double(tmp_r[i][j]) + 0.587 * double(tmp_g[i][j]) + 0.114 * double(tmp_b[i][j]);
            tmp_Cb[i][j] = -0.1687 * double(tmp_r[i][j]) - 0.3313 * double(tmp_g[i][j]) + 0.5 * double(tmp_b[i][j]) + 128;
            tmp_Cr[i][j] = 0.5 * double(tmp_r[i][j]) - 0.4187 * double(tmp_g[i][j]) - 0.0813 * double(tmp_b[i][j]) + 128;
        }
    }

    vector<vector<int>> w0(p_size, vector<int>(p_size, 1));
    for (int i = 0; i < p_size; ++i)
    {
        for (int j = 0; j < p_size; ++j)
        {
            if (fmod(i + j, 2) == 0)
            {
                w0[i][j] = -1;
            }
        }
    }
    cout << "---------w0---------" << p_size << "\n";

    vector<vector<int>> data(32, vector<int>(32, 0));
    vector<int> data2(1024, 0);
    int data_ptr = 0;
    for (int i = 0; i < 32; ++i)
    {
        for (int j = 0; j < 32; ++j)
        {
            srand(time(nullptr));
            data[i][j] = rand() % 2;
            data2[data_ptr] = data[i][j];
            ++data_ptr;
        }
    }
    cout << "---------data---------" << 32 << "\n";

    // 38*38=1444=1023/11==93*15(1395)+1+ 48
    vector<vector<int>> data_hamming(b_size, vector<int>(b_size, 0));
    vector<int> data_temp1(1396, 0);
    vector<int> data_temp2(15, 0);
    vector<int> data_temp(1444, 0);
    data_ptr = 0;
    for (int i = 0; i < 93; ++i)
    {
        int data_j = 11 * i;
        data_temp2[2] = data2[data_j];
        data_temp2[4] = data2[data_j + 1];
        data_temp2[5] = data2[data_j + 2];
        data_temp2[6] = data2[data_j + 3];
        data_temp2[8] = data2[data_j + 4];
        data_temp2[9] = data2[data_j + 5];
        data_temp2[10] = data2[data_j + 6];
        data_temp2[11] = data2[data_j + 7];
        data_temp2[12] = data2[data_j + 8];
        data_temp2[13] = data2[data_j + 9];
        data_temp2[14] = data2[data_j + 10];
        data_temp2[0] = (data_temp2[2] + data_temp2[4] + data_temp2[6] + data_temp2[8] + data_temp2[10] + data_temp2[12] + data_temp2[14]) % 2;
        data_temp2[1] = (data_temp2[2] + data_temp2[5] + data_temp2[6] + data_temp2[9] + data_temp2[10] + data_temp2[13] + data_temp2[14]) % 2;
        data_temp2[3] = (data_temp2[4] + data_temp2[5] + data_temp2[6] + data_temp2[12] + data_temp2[13] + data_temp2[14]) % 2;
        data_temp2[7] = (data_temp2[8] + data_temp2[9] + data_temp2[10] + data_temp2[11] + data_temp2[12] + data_temp2[13] + data_temp2[14]) % 2;

        data_j = 15 * i;
        for (int j = 0; j < 15; ++j)
        {
            data_temp1[data_j + j] = data_temp2[j];
        }
    }
    data_temp1[1395] = data2[1023];

    for (int i = 0; i < 48; ++i)
    {
        data_temp[i] = 0;
    }

    for (int j = 48; j < 1444; ++j)
    {
        data_temp[j] = data_temp1[j - 48];
    }

    for (int i = 0; i < 38; ++i)
    {
        for (int j = 0; j < 38; ++j)
        {
            data_hamming[i][j] = data_temp[i * 38 + j];
        }
    }

    vector<vector<int>> w(w_size, vector<int>(w_size, 0));
    for (int i = 0; i < b_size; ++i)
    {
        for (int j = 0; j < b_size; ++j)
        {
            if (data_hamming[i][j] == 0)
            {
                for (int k = 0; k < p_size; ++k)
                {
                    for (int kk = 0; kk < p_size; ++kk)
                    {
                        w[i * p_size + k][j * p_size + kk] = w0[k][kk];
                    }
                }
            }
            else
            {
                for (int k = 0; k < p_size; ++k)
                {
                    for (int kk = 0; kk < p_size; ++kk)
                    {
                        w[i * p_size + k][j * p_size + kk] = -w0[k][kk];
                    }
                }
            }
        }
    }
    cout << "---------w---------" << w_size << "\n";

    vector<vector<int>> K(w_size, vector<int>(w_size, 0));
    for (int i = 0; i < w_size; ++i)
    {
        for (int j = 0; j < w_size; ++j)
        {
            K[i][j] = (rand() % 2) * 2 - 1;
        }
    }
    cout << "---------K---------" << w_size << "\n";

    vector<vector<int>> w_K_0(w_size, vector<int>(w_size, 0));
    for (int i = 0; i < w_size; ++i)
    {
        for (int j = 0; j < w_size; ++j)
        {
            w_K_0[i][j] = w[i][j] * K[i][j];
        }
    }

    vector<vector<int>> w_K(w_size * 2, vector<int>(w_size * 2, 0));
    for (int i = 0; i < w_size * 2; ++i)
    {
        for (int j = 0; j < w_size * 2; ++j)
        {
            w_K[i][j] = w_K_0[i / 2][j / 2];
        }
    }

    vector<vector<int>> wUD(w_size * 2, vector<int>(w_size * 2, 0));
    for (int i = 0; i < w_size * 2; ++i)
    {
        for (int j = 0; j < w_size * 2; ++j)
        {
            wUD[i][j] = w_K[w_size * 2 - 1 - i][j];
        }
    }

    vector<vector<int>> wLR(w_size * 2, vector<int>(w_size * 2, 0));
    for (int i = 0; i < w_size * 2; ++i)
    {
        for (int j = 0; j < w_size * 2; ++j)
        {
            wLR[i][j] = w_K[i][w_size * 2 - 1 - j];
        }
    }

    vector<vector<int>> wUDLR(w_size * 2, vector<int>(w_size * 2, 0));
    for (int i = 0; i < w_size * 2; ++i)
    {
        for (int j = 0; j < w_size * 2; ++j)
        {
            wUDLR[i][j] = w_K[w_size * 2 - 1 - i][w_size * 2 - 1 - j];
        }
    }

    vector<vector<int>> macro(macro_size, vector<int>(macro_size, 0));
    for (int i = 0; i < w_size * 2; ++i)
    {
        for (int j = 0; j < w_size * 2; ++j)
        {
            macro[i][j] = w_K[i][j];
        }
    }
    for (int i = 0; i < w_size * 2; ++i)
    {
        for (int j = 0; j < w_size * 2; ++j)
        {
            macro[i + w_size * 2][j] = wUD[i][j];
        }
    }
    for (int i = 0; i < w_size * 2; ++i)
    {
        for (int j = 0; j < w_size * 2; ++j)
        {
            macro[i][j + w_size * 2] = wLR[i][j];
        }
    }
    for (int i = 0; i < w_size * 2; ++i)
    {
        for (int j = 0; j < w_size * 2; ++j)
        {
            macro[i + w_size * 2][j + w_size * 2] = wUDLR[i][j];
        }
    }
    cout << "---------macro---------" << w_size * 2 * 2 << "\n";

    int temp_left = len - macro_size * Lnum;

    vector<vector<int>> W(len, vector<int>(len, 0));
    for (int i = 0; i < macro_size * Lnum; ++i)
    {
        for (int j = 0; j < macro_size * Wnum; ++j)
        {
            W[i][j] = macro[fmod(i, macro_size)][fmod(j, macro_size)];
        }
    }
    cout << "---------W(a part of it)---------" << macro_size * Lnum << "\n";

    for (int i = 0; i < temp_left; ++i)
    {
        for (int j = 0; j < macro_size * Lnum; ++j)
        {
            W[j][i + macro_size * Lnum] = macro[j % macro_size][i];
        }
    }

    for (int i = 0; i < temp_left; ++i)
    {
        for (int j = 0; j < temp_left; ++j)
        {
            W[i + macro_size * Lnum][j + macro_size * Lnum] = macro[i][j];
        }
    }

    for (int i = 0; i < temp_left; ++i)
    {
        for (int j = 0; j < macro_size * Lnum; ++j)
        {
            W[i + macro_size * Lnum][j] = macro[i][j % macro_size];
        }
    }

    int lamda = alpha;
    vector<vector<double>> Iy(tmp_height, vector<double>(tmp_width, 0));
    vector<vector<double>> Iyw(tmp_height, vector<double>(tmp_width, 0));
    for (int i = 0; i < tmp_height; ++i)
    {
        for (int j = 0; j < tmp_width; ++j)
        {
            Iy[i][j] = tmp_Y[i][j];
            Iyw[i][j] = tmp_Y[i][j] + lamda * W[i][j];
        }
    }

    vector<vector<int>> tmp_cal_r(tmp_height, vector<int>(tmp_width, 0));
    vector<vector<int>> tmp_cal_g(tmp_height, vector<int>(tmp_width, 0));
    vector<vector<int>> tmp_cal_b(tmp_height, vector<int>(tmp_width, 0));

    Mat tmp_mat_dst;
    tmp_mat_dst.create(I0.size(), I0.type());

    for (int i = 0; i < tmp_height; ++i)
    {
        for (int j = 0; j < tmp_width; ++j)
        {
            tmp_cal_r[i][j] = (int)(1 * Iyw[i][j] + 0 * (tmp_Cb[i][j] - 128) + 1.4 * (tmp_Cr[i][j] - 128));
            tmp_cal_g[i][j] = (int)(1 * Iyw[i][j] - 0.343 * (tmp_Cb[i][j] - 128) - 0.711 * (tmp_Cr[i][j] - 128));
            tmp_cal_b[i][j] = (int)(1 * Iyw[i][j] + 1.765 * (tmp_Cb[i][j] - 128) + 0 * (tmp_Cr[i][j] - 128));
            tmp_mat_dst.at<Vec3b>(i, j)[0] = tmp_cal_b[i][j];
            tmp_mat_dst.at<Vec3b>(i, j)[1] = tmp_cal_g[i][j];
            tmp_mat_dst.at<Vec3b>(i, j)[2] = tmp_cal_r[i][j];
        }
    }

    imwrite("./output/example.png", tmp_mat_dst);

    int d = 3;

    Mat I2 = imread("./output/example.png");
    cout << "---------loading---------\n";
    if (I2.empty() || !I2.data)
    {
        cout << "Loading output image failed!" << endl;
        return -1;
    }
    else
    {
        cout << "Loading output image succsessfully." << endl;
    }
    tmp_height = I2.rows;
    tmp_width = I2.rows;

    for (int i = 0; i < tmp_height; ++i)
    {
        for (int j = 0; j < tmp_width; ++j)
        {
            tmp_b[i][j] = I2.at<Vec3b>(i, j)[0];
            tmp_g[i][j] = I2.at<Vec3b>(i, j)[1];
            tmp_r[i][j] = I2.at<Vec3b>(i, j)[2];
        }
    }
    for (int i = 0; i < tmp_height; ++i)
    {
        for (int j = 0; j < tmp_width; ++j)
        {
            tmp_Y[i][j] = 0.299 * double(tmp_r[i][j]) + 0.587 * double(tmp_g[i][j]) + 0.114 * double(tmp_b[i][j]);
            tmp_Cb[i][j] = -0.1687 * double(tmp_r[i][j]) - 0.3313 * double(tmp_g[i][j]) + 0.5 * double(tmp_b[i][j]) + 128;
            tmp_Cr[i][j] = 0.5 * double(tmp_r[i][j]) - 0.4187 * double(tmp_g[i][j]) - 0.0813 * double(tmp_b[i][j]) + 128;
        }
    }

    // cropping
    double temp_cut_rate = 0.5;
    int temp_cut_len = int(2048 * sqrt(1 - temp_cut_rate));
    cout << "------------------\n";
    cout << "---------temp_cut_rate---------" << temp_cut_rate << "\n";
    cout << "---------temp_cut_len---------" << temp_cut_len << "\n";
    cout << "------------------\n";
    tmp_height = temp_cut_len;
    tmp_width = temp_cut_len;

    vector<vector<double>> tmp_Y_cut(tmp_height, vector<double>(tmp_width, 0));
    vector<int> temp_rand_cut_start_X(10, 0);
    vector<int> temp_rand_cut_start_Y(10, 0);
    for (int i = 0; i < 10; ++i)
    {
        int random_number = rand() % (2048 - tmp_height + 1);
        temp_rand_cut_start_X[i] = random_number;
    }
    for (int i = 0; i < 10; ++i)
    {
        int random_number = rand() % (2048 - tmp_height + 1);
        temp_rand_cut_start_Y[i] = random_number;
    }

    for (int tmp_round_i = 0; tmp_round_i < 1; ++tmp_round_i)
    {
        for (int i = 0; i < tmp_height; ++i)
        {
            for (int j = 0; j < tmp_width; ++j)
            {
                tmp_Y_cut[i][j] = tmp_Y[i + temp_rand_cut_start_X[tmp_round_i]][j + temp_rand_cut_start_Y[tmp_round_i]];
            }
        }

        cout << "temp_rand_cut_start_X: " << temp_rand_cut_start_X[tmp_round_i] << "  temp_rand_cut_start_Y: " << temp_rand_cut_start_Y[tmp_round_i] << "\n";

        vector<vector<double>> localMean(tmp_height, vector<double>(tmp_width, 0));
        for (int i = -d / 2; i < tmp_height - d / 2; i++)
        {
            for (int j = -d / 2; j < tmp_width - d / 2; j++)
            {
                double filtertotal = 0;
                double filternumber = 0;
                for (int k = 0; k < d; k++)
                {
                    for (int l = 0; l < d; l++)
                    {
                        if (i + k >= 0 && i + k < tmp_height && j + l >= 0 && j + l < tmp_width)
                        {
                            filtertotal += tmp_Y_cut[i + k][j + l];
                            filternumber += 1.0;
                        }
                    }
                }
                localMean[i + (d / 2)][j + (d / 2)] = filtertotal / filternumber;
            }
        }

        vector<vector<double>> tmp_Y_cut_square(tmp_height, vector<double>(tmp_width, 0));
        for (int i = 0; i < tmp_height; ++i)
        {
            for (int j = 0; j < tmp_width; ++j)
            {
                tmp_Y_cut_square[i][j] = tmp_Y_cut[i][j] * tmp_Y_cut[i][j];
            }
        }

        vector<vector<double>> localMean_square(tmp_height, vector<double>(tmp_width, 0));
        for (int i = 0; i < tmp_height; ++i)
        {
            for (int j = 0; j < tmp_width; ++j)
            {
                localMean_square[i][j] = localMean[i][j] * localMean[i][j];
            }
        }

        vector<vector<double>> tmp_localVar0(tmp_height, vector<double>(tmp_width, 0));
        for (int i = -d / 2; i < tmp_height - d / 2; i++)
        {
            for (int j = -d / 2; j < tmp_width - d / 2; j++)
            {
                double filtertotal = 0;
                double filternumber = 0;
                for (int k = 0; k < d; k++)
                {
                    for (int l = 0; l < d; l++)
                    {
                        if (i + k >= 0 && i + k < tmp_height && j + l >= 0 && j + l < tmp_width)
                        {
                            filtertotal += tmp_Y_cut_square[i + k][j + l];
                            filternumber += 1.0;
                        }
                    }
                }
                tmp_localVar0[i + (d / 2)][j + (d / 2)] = filtertotal / filternumber;
            }
        }

        vector<vector<double>> localVar(tmp_height, vector<double>(tmp_width, 0));
        for (int i = 0; i < tmp_height; ++i)
        {
            for (int j = 0; j < tmp_width; ++j)
            {
                localVar[i][j] = tmp_localVar0[i][j] - localMean_square[i][j];
            }
        }

        double noi = 0;
        for (int i = 0; i < tmp_height; ++i)
        {
            for (int j = 0; j < tmp_width; ++j)
            {
                noi += localVar[i][j];
            }
        }
        noi = noi / tmp_height / tmp_width;
        // cout << "---------noi---------" << noi << "\n";

        for (int i = 0; i < tmp_height; ++i)
        {
            for (int j = 0; j < tmp_width; ++j)
            {
                if (localVar[i][j] < noi)
                {
                    localVar[i][j] = noi;
                }
            }
        }
        // cout << "---------localVar---------" << tmp_height << "\n";

        vector<vector<double>> w_est(tmp_height, vector<double>(tmp_width, 0));
        for (int i = 0; i < tmp_height; ++i)
        {
            for (int j = 0; j < tmp_width; ++j)
            {
                w_est[i][j] = (tmp_Y_cut[i][j] - localMean[i][j]) * noi / localVar[i][j];
            }
        }
        cout << "---------w_est---------" << tmp_height << "\n";

        vector<double> temp_row_sym(500, 0);

        for (int i = tmp_height / 2 - 250; i < tmp_height / 2 + 250; ++i)
        {
            for (int j = 2; j < 600; ++j)
            {
                for (int k = 0; k < 100; ++k)
                {
                    temp_row_sym[i - (tmp_height / 2 - 250)] += w_est[i - k + 1][j] * w_est[i + k][j];
                }
            }
        }

        double temp_row_sym_sum = 0;
        for (int i = 0; i < 500; ++i)
        {
            temp_row_sym_sum += temp_row_sym[i];
        }
        temp_row_sym_sum /= 100;

        int temp_row_sym_count = 0;
        vector<double> temp_row_sym_res(500, 0);
        vector<double> temp_row_sym_num(500, 0);

        for (int i = 0; i < 500; ++i)
        {
            if (temp_row_sym[i] > temp_row_sym_sum)
            {
                temp_row_sym_res[temp_row_sym_count] = temp_row_sym[i];
                temp_row_sym_num[temp_row_sym_count] = i;
                ++temp_row_sym_count;
            }
        }

        int temp_row_max = 0;
        int temp_row_max_num = 0;
        for (int i = 0; i < temp_row_sym_count; ++i)
        {
            if (temp_row_sym_res[i] > temp_row_max)
            {
                temp_row_max = temp_row_sym_res[i];
                temp_row_max_num = temp_row_sym_num[i];
            }
        }
        int temp_row_flag0 = 0;
        int temp_row_flag1 = 0;
        int temp_row_flag0_value = 0;
        int temp_row_flag1_value = 0;
        for (int i = 0; i < temp_row_sym_count; ++i)
        {
            if (temp_row_sym_num[i] == temp_row_max_num - 1)
            {
                temp_row_flag0 = 1;
                temp_row_flag0_value = temp_row_sym_res[i];
            }
            if (temp_row_sym_num[i] == temp_row_max_num + 1)
            {
                temp_row_flag1 = 1;
                temp_row_flag1_value = temp_row_sym_res[i];
            }
        }
        int temp_row_real = temp_row_max_num + tmp_height / 2 - 250;
        cout << "sym_row: " << temp_row_max_num << "\n";
        cout << "temp_row_real: " << temp_row_real % 152 << "\n";

        vector<double> temp_col_sym(500, 0);

        for (int i = tmp_height / 2 - 250; i < tmp_height / 2 + 250; ++i)
        {
            for (int j = 0; j < 600; ++j)
            {
                for (int k = 2; k < 100; ++k)
                {
                    temp_col_sym[i - (tmp_height / 2 - 250)] += w_est[j][i - k + 1] * w_est[j][i + k];
                }
            }
        }

        double temp_col_sym_sum = 0;
        for (int i = 0; i < 500; ++i)
        {
            temp_col_sym_sum += temp_col_sym[i];
        }
        temp_col_sym_sum /= 100;

        int temp_col_sym_count = 0;
        vector<double> temp_col_sym_res(500, 0);
        vector<double> temp_col_sym_num(500, 0);

        for (int i = 0; i < 500; ++i)
        {
            if (temp_col_sym[i] > temp_col_sym_sum)
            {
                temp_col_sym_res[temp_col_sym_count] = temp_col_sym[i];
                temp_col_sym_num[temp_col_sym_count] = i;
                ++temp_col_sym_count;
            }
        }

        int temp_col_max = 0;
        int temp_col_max_num = 0;
        for (int i = 0; i < temp_col_sym_count; ++i)
        {
            if (temp_col_sym_res[i] > temp_col_max)
            {
                temp_col_max = temp_col_sym_res[i];
                temp_col_max_num = temp_col_sym_num[i];
            }
        }
        int temp_col_flag0 = 0;
        int temp_col_flag1 = 0;
        int temp_col_flag0_value = 0;
        int temp_col_flag1_value = 0;
        for (int i = 0; i < temp_col_sym_count; ++i)
        {
            if (temp_col_sym_num[i] == temp_col_max_num - 1)
            {
                temp_col_flag0 = 1;
                temp_col_flag0_value = temp_col_sym_res[i];
            }
            if (temp_col_sym_num[i] == temp_col_max_num + 1)
            {
                temp_col_flag1 = 1;
                temp_col_flag1_value = temp_col_sym_res[i];
            }
        }
        int temp_col_real = temp_col_max_num + tmp_height / 2 - 250;
        cout << "sym_col: " << temp_col_max_num << "\n";
        cout << "temp_col_real: " << temp_col_real % 152 << "\n";

        int temp_row_start = temp_row_real % 152 + 1;
        int temp_col_start = temp_col_real % 152 + 1;

        vector<vector<double>> temp_area_total(304, vector<double>(304, 0));
        for (int i = 0; i < 304; ++i)
        {
            for (int j = 0; j < 304; ++j)
            {
                temp_area_total[i][j] = w_est[temp_row_start + i][temp_col_start + j];
            }
        }

        vector<vector<double>> temp_area1(76, vector<double>(76, 0));
        for (int i = 0; i < 76; ++i)
        {
            for (int j = 0; j < 76; ++j)
            {
                temp_area1[i][j] = temp_area_total[i * 2][j * 2] + temp_area_total[i * 2 + 1][j * 2] + temp_area_total[i * 2][j * 2 + 1] + temp_area_total[i * 2 + 1][j * 2 + 1];
                temp_area1[i][j] = temp_area1[i][j] * K[i][j];
            }
        }
        int temp_area1_local_true = 0;
        for (int i = 0; i < 38; ++i)
        {
            for (int j = 0; j < 38; ++j)
            {
                if ((temp_area1[2 * i][2 * j] > temp_area1[2 * i + 1][2 * j] && temp_area1[2 * i][2 * j] > temp_area1[2 * i][2 * j + 1] && temp_area1[2 * i + 1][2 * j + 1] > temp_area1[2 * i + 1][2 * j] && temp_area1[2 * i + 1][2 * j + 1] > temp_area1[2 * i][2 * j + 1]) || (temp_area1[2 * i][2 * j] < temp_area1[2 * i + 1][2 * j] && temp_area1[2 * i][2 * j] < temp_area1[2 * i][2 * j + 1] && temp_area1[2 * i + 1][2 * j + 1] < temp_area1[2 * i + 1][2 * j] && temp_area1[2 * i + 1][2 * j + 1] < temp_area1[2 * i][2 * j + 1]))
                {
                    temp_area1_local_true++;
                }
            }
        }
        vector<vector<double>> temp_area2(76, vector<double>(76, 0));
        for (int i = 0; i < 76; ++i)
        {
            for (int j = 0; j < 76; ++j)
            {
                temp_area2[i][j] = temp_area_total[(i + 76) * 2][j * 2] + temp_area_total[(i + 76) * 2 + 1][j * 2] + temp_area_total[(i + 76) * 2][j * 2 + 1] + temp_area_total[(i + 76) * 2 + 1][j * 2 + 1];
                temp_area2[i][j] = temp_area2[i][j] * K[i][j];
            }
        }
        int temp_area2_local_true = 0;
        for (int i = 0; i < 38; ++i)
        {
            for (int j = 0; j < 38; ++j)
            {
                if ((temp_area2[2 * i][2 * j] > temp_area2[2 * i + 1][2 * j] && temp_area2[2 * i][2 * j] > temp_area2[2 * i][2 * j + 1] && temp_area2[2 * i + 1][2 * j + 1] > temp_area2[2 * i + 1][2 * j] && temp_area2[2 * i + 1][2 * j + 1] > temp_area2[2 * i][2 * j + 1]) || (temp_area2[2 * i][2 * j] < temp_area2[2 * i + 1][2 * j] && temp_area2[2 * i][2 * j] < temp_area2[2 * i][2 * j + 1] && temp_area2[2 * i + 1][2 * j + 1] < temp_area2[2 * i + 1][2 * j] && temp_area2[2 * i + 1][2 * j + 1] < temp_area2[2 * i][2 * j + 1]))
                {
                    temp_area2_local_true++;
                }
            }
        }
        vector<vector<double>> temp_area3(76, vector<double>(76, 0));
        for (int i = 0; i < 76; ++i)
        {
            for (int j = 0; j < 76; ++j)
            {
                temp_area3[i][j] = temp_area_total[(i + 76) * 2][(j + 76) * 2] + temp_area_total[(i + 76) * 2 + 1][(j + 76) * 2] + temp_area_total[(i + 76) * 2][(j + 76) * 2 + 1] + temp_area_total[(i + 76) * 2 + 1][(j + 76) * 2 + 1];
                temp_area3[i][j] = temp_area3[i][j] * K[i][j];
            }
        }
        int temp_area3_local_true = 0;
        for (int i = 0; i < 38; ++i)
        {
            for (int j = 0; j < 38; ++j)
            {
                if ((temp_area3[2 * i][2 * j] > temp_area3[2 * i + 1][2 * j] && temp_area3[2 * i][2 * j] > temp_area3[2 * i][2 * j + 1] && temp_area3[2 * i + 1][2 * j + 1] > temp_area3[2 * i + 1][2 * j] && temp_area3[2 * i + 1][2 * j + 1] > temp_area3[2 * i][2 * j + 1]) || (temp_area3[2 * i][2 * j] < temp_area3[2 * i + 1][2 * j] && temp_area3[2 * i][2 * j] < temp_area3[2 * i][2 * j + 1] && temp_area3[2 * i + 1][2 * j + 1] < temp_area3[2 * i + 1][2 * j] && temp_area3[2 * i + 1][2 * j + 1] < temp_area3[2 * i][2 * j + 1]))
                {
                    temp_area3_local_true++;
                }
            }
        }
        vector<vector<double>> temp_area4(76, vector<double>(76, 0));
        for (int i = 0; i < 76; ++i)
        {
            for (int j = 0; j < 76; ++j)
            {
                temp_area4[i][j] = temp_area_total[i * 2][(j + 76) * 2] + temp_area_total[i * 2 + 1][(j + 76) * 2] + temp_area_total[i * 2][(j + 76) * 2 + 1] + temp_area_total[i * 2 + 1][(j + 76) * 2 + 1];
                temp_area4[i][j] = temp_area4[i][j] * K[i][j];
            }
        }
        int temp_area4_local_true = 0;
        for (int i = 0; i < 38; ++i)
        {
            for (int j = 0; j < 38; ++j)
            {
                if ((temp_area4[2 * i][2 * j] > temp_area4[2 * i + 1][2 * j] && temp_area4[2 * i][2 * j] > temp_area4[2 * i][2 * j + 1] && temp_area4[2 * i + 1][2 * j + 1] > temp_area4[2 * i + 1][2 * j] && temp_area4[2 * i + 1][2 * j + 1] > temp_area4[2 * i][2 * j + 1]) || (temp_area4[2 * i][2 * j] < temp_area4[2 * i + 1][2 * j] && temp_area4[2 * i][2 * j] < temp_area4[2 * i][2 * j + 1] && temp_area4[2 * i + 1][2 * j + 1] < temp_area4[2 * i + 1][2 * j] && temp_area4[2 * i + 1][2 * j + 1] < temp_area4[2 * i][2 * j + 1]))
                {
                    temp_area4_local_true++;
                }
            }
        }

        vector<vector<double>> temp_data_est(38, vector<double>(38, 0));

        if (temp_area1_local_true > temp_area2_local_true && temp_area1_local_true > temp_area3_local_true && temp_area1_local_true > temp_area4_local_true)
        {
            for (int i = 0; i < 38; ++i)
            {
                for (int j = 0; j < 38; ++j)
                {
                    temp_data_est[i][j] = temp_area1[i * 2][j * 2] - temp_area1[i * 2 + 1][j * 2] - temp_area1[i * 2][j * 2 + 1] + temp_area1[i * 2 + 1][j * 2 + 1];
                    if (temp_data_est[i][j] >= 0)
                    {
                        temp_data_est[i][j] = 1;
                    }
                    else
                    {
                        temp_data_est[i][j] = 0;
                    }
                }
            }
        }
        else if (temp_area2_local_true > temp_area1_local_true && temp_area2_local_true > temp_area3_local_true && temp_area1_local_true > temp_area4_local_true)
        {
            for (int i = 0; i < 38; ++i)
            {
                for (int j = 0; j < 38; ++j)
                {
                    temp_data_est[i][j] = temp_area2[i * 2][j * 2] - temp_area2[i * 2 + 1][j * 2] - temp_area2[i * 2][j * 2 + 1] + temp_area2[i * 2 + 1][j * 2 + 1];
                    if (temp_data_est[i][j] >= 0)
                    {
                        temp_data_est[i][j] = 1;
                    }
                    else
                    {
                        temp_data_est[i][j] = 0;
                    }
                }
            }
        }
        else if (temp_area3_local_true > temp_area1_local_true && temp_area3_local_true > temp_area2_local_true && temp_area1_local_true > temp_area4_local_true)
        {
            for (int i = 0; i < 38; ++i)
            {
                for (int j = 0; j < 38; ++j)
                {
                    temp_data_est[i][j] = temp_area3[i * 2][j * 2] - temp_area3[i * 2 + 1][j * 2] - temp_area3[i * 2][j * 2 + 1] + temp_area3[i * 2 + 1][j * 2 + 1];
                    if (temp_data_est[i][j] >= 0)
                    {
                        temp_data_est[i][j] = 1;
                    }
                    else
                    {
                        temp_data_est[i][j] = 0;
                    }
                }
            }
        }
        else if (temp_area4_local_true > temp_area1_local_true && temp_area4_local_true > temp_area2_local_true && temp_area4_local_true > temp_area3_local_true)
        {
            for (int i = 0; i < 38; ++i)
            {
                for (int j = 0; j < 38; ++j)
                {
                    temp_data_est[i][j] = temp_area4[i * 2][j * 2] - temp_area4[i * 2 + 1][j * 2] - temp_area4[i * 2][j * 2 + 1] + temp_area4[i * 2 + 1][j * 2 + 1];
                    if (temp_data_est[i][j] >= 0)
                    {
                        temp_data_est[i][j] = 1;
                    }
                    else
                    {
                        temp_data_est[i][j] = 0;
                    }
                }
            }
        }
    }
    return 0;
}