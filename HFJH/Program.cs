using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Data;
using NSMatrix;
using System.Data.OleDb;

namespace HFJH
{
    class Program
    {
        static void Main(string[] args)
        {
            double[] pic_coor_x = new double[132];//像平面x坐标
            double[] pic_coor_y = new double[132];//像平面y坐标
            double[] con_coor_x = new double[132];//物方控制点坐标
            double[] con_coor_y = new double[132];
            double[] con_coor_z = new double[132];
            double[] x_coor = new double[135];//像点x坐标
            double[] y_coor = new double[132];//像点y坐标
            double[] x_0 = new double[132];//x，y近似值
            double[] y_0 = new double[132];

            double f = 4547.9352;//主距
            double xs, ys, zs;//投影中心的地面坐标系坐标
            double p = 0, o = 0, k = 0;//估计欧拉角
            double Xs0 = 3370,Ys0 = -140,Zs0 = 90;//估计投影中心坐标初值

          
            double a1, a2, a3, b1, b2, b3, c1, c2, c3;//方向余弦
            
            string picx = @"D:\DATA\像点x.txt";
            string picy = @"D:\DATA\像点y.txt";
            string conx = @"D:\DATA\控制点x.txt";
            string cony = @"D:\DATA\控制点y.txt";
            string conz = @"D:\DATA\控制点z.txt";
            Text(picx, x_coor);
            Text(picy, y_coor);
            Text(conx,con_coor_x);
            Text(cony,con_coor_y);
            Text(conz,con_coor_z);
            Console.WriteLine("数据读取成功！");
            Console.WriteLine("按任意键继续！");
            Console.ReadKey();
            Transpic(x_coor, pic_coor_x, y_coor, pic_coor_y);
            
            double[] R = new double[9];//旋转矩阵
            a1 =  Math.Cos(p) * Math.Cos(k) - Math.Sin(p) * Math.Sin(o) * Math.Sin(k);
            a2 = -Math.Cos(p) * Math.Sin(k) - Math.Sin(p) * Math.Sin(o) * Math.Cos(k);
            a3 = -Math.Sin(p) * Math.Cos(o);
            b1 = Math.Cos(o) * Math.Sin(k);
            b2 = Math.Cos(o) * Math.Cos(k);
            b3 = -Math.Sin(o);
            c1 = Math.Sin(p) * Math.Cos(k) + Math.Cos(p) * Math.Sin(o) * Math.Sin(k);
            c2 = -Math.Sin(p) * Math.Sin(k) + Math.Cos(p) * Math.Sin(o) * Math.Cos(k);
            c3 = Math.Cos(p) * Math.Cos(o);
            R[0] = a1;//将方向余弦存入变换矩阵
            R[1] = a2;
            R[2] = a3;
            R[3] = b1;
            R[4] = b2;
            R[5] = b3;
            R[6] = c1;
            R[7] = c2;
            R[8] = c3;
            for(int i =0;i<132;i++)//求得x，y的近似坐标
            {
                x_0[i] = (-f) * (a1 * (con_coor_x[i] - Xs0) + b1 * (con_coor_y[i] - Ys0) + c1 * (con_coor_z[i] - Zs0)) / (a3 * (con_coor_x[i] - Xs0) + b3 * (con_coor_y[i] - Ys0) + c3 * (con_coor_z[i] - Zs0));
                y_0[i]=  (-f) * (a2 * (con_coor_x[i] - Xs0) + b2 * (con_coor_y[i] - Ys0) + c2 * (con_coor_z[i] - Zs0)) / (a3 * (con_coor_x[i] - Xs0) + b3 * (con_coor_y[i] - Ys0) + c3 * (con_coor_z[i] - Zs0));
            }
            double[,] b = new double[264, 6];//存放系数矩阵的值
            double[] _Z = new double[132];//系数计算中的z
            for (int i = 0; i < 132; i++)
            {
                _Z[i] = (a3 * (con_coor_x[i] - Xs0) + b3 * (con_coor_y[i] - Ys0) + c3 * (con_coor_z[i] - Zs0));
            }
            for (int i = 0; i < 132; i++)

            {

                //计算系数矩阵

                b[2 * i, 0] = (a1 * f + a2 * pic_coor_x[i]) / _Z[i];

                b[2 * i, 1] = (b1 * f + b3 * pic_coor_x[i]) / _Z[i];

                b[2 * i, 2] = (c1 * f + c3 * pic_coor_x[i]) / _Z[i];

                b[2 * i, 3] = pic_coor_y[i] * Math.Sin(o) - ((pic_coor_x[i] / f) * (pic_coor_x[i] * Math.Cos(k) - pic_coor_y[i] * Math.Sin(k)) + f * Math.Cos(k)) * Math.Cos(o);

                b[2 * i, 4] = -f * Math.Sin(k) - (pic_coor_x[i] / f) * (pic_coor_x[i] * Math.Sin(k) + pic_coor_y[i] * Math.Cos(k));

                b[2 * i, 5] = pic_coor_y[i];



                b[2 * i + 1, 0] = (a2 * f + a3 * pic_coor_y[i]) / _Z[i];

                b[2 * i + 1, 1] = (b2 * f +b3 * pic_coor_y[i]) / _Z[i];

                b[2 * i + 1, 2] = (c2 * f + c3 * pic_coor_y[i]) / _Z[i];

                b[2 * i + 1, 3] = -pic_coor_x[i] * Math.Sin(p) - ((pic_coor_x[i] / f) * (pic_coor_x[i] * Math.Cos(k) - pic_coor_y[i] * Math.Sin(k)) - f * Math.Sin(k)) * Math.Cos(o);

                b[2 * i + 1, 4] = -f * Math.Cos(k) - (pic_coor_y[i] / f) * (pic_coor_x[i] * Math.Sin(k) + pic_coor_y[i] * Math.Cos(k));

                b[2 * i + 1, 5] = -pic_coor_x[i];
            }//系数矩阵元素计算
            Matrix B = new Matrix(b);//系数矩阵B
            double[,] l = new double[264, 1];
            for (int i = 0; i < 132; i++)//计算常数L
            {
                l[2 * i, 0] = pic_coor_x[i] - x_0[i];

                l[2 * i + 1, 0] = pic_coor_y[i] - y_0[i];

            }
            Matrix L = new Matrix(l);//常数矩阵L
            Matrix B_t = B.Transpose();
            Matrix C = (B_t * B);//系数矩阵的转置与系数矩阵相乘
            C.InvertGaussJordan();
          
            Matrix D = C * B_t;
            Matrix E = D * L;
            double[] res = new double[6];
            E.Show_Mat(res);//E为最终的改正数
            double Xs = Xs0 + res[0];
            double Ys = Ys0 + res[1];
            double Zs = Zs0 + res[2];
            double phi = p + res[3];
            double omiga = o + res[4];
            double kappa = k + res[5];
            Console.Write("Xs改正后为：");
            Console.WriteLine(Xs);
            Console.Write("Ys改正后为：");
            Console.WriteLine(Ys);
            Console.Write("Zs改正后为：");
            Console.WriteLine(Zs);
            Console.Write("phi改正后为：");
            Console.WriteLine(phi);
            Console.Write("omiga改正后为：");
            Console.WriteLine(omiga);
            Console.Write("kappa改正后为：");
            Console.WriteLine(kappa);
            Console.ReadKey();


        }
        static void Transpic(double[] a,double[] b,double[] c,double[] d)//像点坐标转换像平面坐标
        {
            int h = 4008;//像片高
            int w = 5344;//像片宽
            double x0 = 47.4857;//像主点相对影像中心的位置
            double y0 = 12.0276;
            for (int i=0; i<132; i++)//像点坐标转化为像平面坐标
            {
                b[i] = (a[i] -( w * 0.5)) - x0;//a是像点x坐标，b是像平面x坐标
                d[i] = (h * 0.5) - c[i] - y0;//c是像点y坐标，d是像平面y坐标
            }
        }
        static void Text(string src,double[] a)//文件读取并存入数组
        {
            StreamReader streamReader = new StreamReader(src);
            String[] fileLines = System.IO.File.ReadAllLines(src);
            for (int i = 0; i < 132; i++)
            {
                a[i] = Convert.ToDouble(fileLines[i]);
            }
        }
    }
}
