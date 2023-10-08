using System;

namespace Com_Methods
{
    /// <summary>
    /// класс, реализующий решатель на базе LU-разложения
    /// </summary>
    class LU_Decomposition
    {
        //матрица-хранилище для LU-разложения
        public Matrix LU { set; get; }
        //порядок следования строк
        int[] P { set; get; }

        //конструктор по умолчанию
        public LU_Decomposition(){}

        //реализация LU-разложения
        public LU_Decomposition(Matrix A)
        {
            //хранилище для верхней и нижней треугольных матриц
            LU = new Matrix(A.M, A.N);
            //копирование исходной матрицы
            LU.Copy(A);

            //инициализация матрицы перестановок строк
            P = new int[A.M];
            for (int i = 0; i < P.Length; i++) P[i] = i;

            //построение верхней треугольной матрицы
            for (int i = 0; i < LU.M - 1; i++)
            {
                //находим ведущий элемент в i-том столбце
                int I = Gauss_Method.Find_Main_Element(LU, i);

                //если это не диагональ
                if (I != i)
                {
                    //переставляем строки I и i в СЛАУ местами
                    var Row = LU.Elem[I];
                    LU.Elem[I] = LU.Elem[i];
                    LU.Elem[i] = Row;
                    int Index = P[i];
                    P[i] = P[I];
                    P[I] = Index;
                }
                
                //для оставшихся строк выполним умножение слева на матрицу преобразований
                for (int j = i + 1; j < LU.M; j++)
                {
                    double help = LU.Elem[j][i] / LU.Elem[i][i];

                    //для уменьшения ошибок вычислений обнуляемые компоненты занулим явно 
                    LU.Elem[j][i] = 0;

                    //вычитаем элементы строки i из строк от i + 1 до A.M 
                    for (int k = i + 1; k < LU.M; k++)
                    {
                        LU.Elem[j][k] -= help * LU.Elem[i][k];
                    }
                }
            }

            //построение нижней треугольной матрицы
            for (int i = 0; i < A.M; i++)
            {
                for (int j = 0; j < i; j++)
                {
                    double sum_LikUkj = 0;
                    for (int k = 0; k < j; k++)
                    {
                        sum_LikUkj += LU.Elem[i][k] * LU.Elem[k][j];
                    }
                    LU.Elem[i][j] = (A.Elem[P[i]][j] - sum_LikUkj) / LU.Elem[j][j];
                }
            }
        }

        /// <summary>
        /// прямой ход (без деления на диагональ, т.к. L_ii = 1)
        /// </summary>
        void Direct_Way(Vector F, Vector RES)
        {
            for (int i = 0; i < LU.M; i++) RES.Elem[i] = F.Elem[P[i]];

            for (int i = 0; i < LU.M; i++)
            {
                for (int j = 0; j < i; j++)
                {
                    RES.Elem[i] -= LU.Elem[i][j] * RES.Elem[j];
                }
            }
        }

        /// <summary>
        /// обратный ход (по строкам)
        /// </summary>
        void Back_Way(Vector F, Vector RES)
        {
            for (int i = 0; i < LU.M; i++) RES.Elem[i] = F.Elem[i];

            for (int i = F.N - 1; i >= 0; i--)
            {
                if (Math.Abs(LU.Elem[i][i]) < CONST.EPS) throw new Exception("Back Row Substitution: division by 0... ");

                for (int j = i + 1; j < F.N; j++)
                {
                    RES.Elem[i] -= LU.Elem[i][j] * RES.Elem[j];
                }

                RES.Elem[i] /= LU.Elem[i][i];
            }
        }

        /// <summary>
        /// решатель СЛАУ с плотной матрицей на базе LU-разложения
        /// LU разложение предполагается построенным
        /// </summary>
        public Vector Start_Solver(Vector F)
        {
            if (LU == null) throw new Exception("Error: LU-matrix is not built. Constructing of LU-matrix is required...");
            var RES = new Vector(F.N);
            //прямой ход
            Direct_Way(F, RES);
            //обратный ход
            Back_Way(RES, RES);
            return RES;
        }
    }
}