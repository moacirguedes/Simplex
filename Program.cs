using System;
using System.Collections.Generic;

namespace Simplex
{
    class Program
    {
        static void Main(string[] args)
        {
            System.Console.WriteLine("Entre com o tipo da funcao.\n 0: Max \n 1: Min");
            int Min = Convert.ToInt32(Console.ReadLine());
            System.Console.WriteLine("\nEntre com o numero de variaveis:");
            int NumVariaveis = Convert.ToInt32(Console.ReadLine());
            System.Console.WriteLine("\nEntre com o numero de restricoes:");
            int NumRestricoes = Convert.ToInt32(Console.ReadLine());

            double[] Funcao = new double[NumVariaveis];

            double[,] Restricoes = new double[NumRestricoes, NumVariaveis];
            int[] TipoRestricoes = new int[NumRestricoes];
            double[] ValorRestricoes = new double[NumRestricoes];
            int[] TipoVariaveis = new int[NumVariaveis];

            double[,] Base = new double[NumRestricoes, NumRestricoes];
            int[] VarBase = new int[NumRestricoes];
            int[] VarNBase;

            double [,] A;
            int tamA;

            Gauss gauss = new Gauss();
            int it = 0;
            double[] Lambda;
            double[] X;
            double[] xB;
            double[,] aux;
            double[,] bT;
            double[] cB;
            double[] cN;
            int[] entraBase = new int[2];
            double[] nK;
            double[] y;
            int[] saiBase = new int[2];


            for (int i = 0; i < NumVariaveis; i++)
            {
                System.Console.WriteLine("Entre com o valor que multiplica x{0} na funcao: ", i+1);
                Funcao[i] = Convert.ToDouble(Console.ReadLine());
            }

            for (int i = 0; i < NumRestricoes; i++)
            {
                for (int j = 0; j < NumVariaveis; j++)
                {
                    System.Console.WriteLine("Entre com o valor que multipliza x{0} na restricao {1}: ", j+1, i+1);
                    Restricoes[i, j] = Convert.ToDouble(Console.ReadLine());
                }
            }

            for (int i = 0; i < NumRestricoes; i++)
            {
                System.Console.WriteLine("Entre com o tipo da restricao {0}. ([0 =] [-1 >=] [1 <=]): ", i+1);
                TipoRestricoes[i] = Convert.ToInt32(Console.ReadLine());
            }

            for (int i = 0; i < NumRestricoes; i++)
            {
                System.Console.WriteLine("Entre com o valor da restricao {0}: ", i+1);
                ValorRestricoes[i] = Convert.ToDouble(Console.ReadLine());
            }

            for (int i = 0; i < NumVariaveis; i++)
            {
                System.Console.WriteLine("Entre com o tipo de x{0}. ([0 livre] [-1 >=] [1 <=]): ", i+1);
                TipoVariaveis[i] = Convert.ToInt32(Console.ReadLine());
            }
            System.Console.WriteLine();

            // ----------------Verificação de necessidade da Fase 1

            bool FaseUm = false;
            int VarArtificial = 0;

            // Se for Max
            if (Min == 0)
            {
                for (int i = 0; i < NumVariaveis; i++)
                {
                    Funcao[i] *= -1;
                }
            }

            // Se tiver b < 0
            for (int i = 0; i < NumRestricoes; i++)
            {
                if (ValorRestricoes[i] < 0)
                {
                    for (int j = 0; j < NumVariaveis; j++)
                    {
                        Restricoes[i, j] *= -1;
                    }
                    TipoRestricoes[i] *= -1;
                    ValorRestricoes[i] *= -1;
                }
            }
          
            // Verifica necessidade de Fase1 e a quantidade de variaveis artificiais
            List<int> RestricoesErradas = new List<int>();
            List<int> RestricoesIgual = new List<int>();
            List<int> RestricoesMaior = new List<int>();
            for (int i = 0; i < NumRestricoes; i++)
            {
                if (TipoRestricoes[i] != 1)
                {
                    FaseUm = true;
                    RestricoesErradas.Add(i);
                    VarArtificial++;
                    if (TipoRestricoes[i] == 0)
                    {
                        RestricoesIgual.Add(i);
                    }
                    else
                    {
                        RestricoesMaior.Add(i);
                    }
                }
            }

            // Cria a matriz A com as Retrições e as novas Variaveis Auxiliares
            A = Restricoes;
            tamA = NumVariaveis;
            double [] vet = new double[NumRestricoes];
            for (int i = 0 ; i < NumRestricoes ; i++)
            {
                if (RestricoesMaior.Contains(i))
                {
                    for (int j = 0 ; j < NumRestricoes ; j++)
                    {
                        if (i == j)
                        {
                            vet[j] = -1;
                        }
                        else
                        {
                            vet[j] = 0;
                        }
                    }
                }
                else
                {
                    if (RestricoesIgual.Contains(i))
                    {
                        for (int j = 0 ; j < NumRestricoes ; j++)
                        {
                            vet[j] = 0;
                        }
                    }
                    else
                    {
                        for (int j = 0 ; j < NumRestricoes ; j++)
                        {
                            if (i == j)
                            {
                                vet[j] = 1;
                            }
                            else
                            {
                                vet[j] = 0;
                            }
                        }
                    }
                }
                A = MatrizAddColum(NumRestricoes, tamA, A, vet);
                tamA++;
            }

            // Adiciona as colunas das Variaveis Artificiais na matriz A
            for (int i = 0; i < RestricoesErradas.Count ; i++)
            {
                for (int j = 0 ; j < NumRestricoes ; j++)
                {
                    if (RestricoesErradas[i] == j)
                    {
                        vet[j] = 1;
                    }
                    else
                    {
                        vet[j] = 0;
                    }
                }
                A = MatrizAddColum(NumRestricoes, tamA, A, vet);
                tamA++;
            }

            // Cria Matriz Base = identidade
            for (int i = 0 ; i < NumRestricoes ; i++)
            {
                for (int j = 0 ; j < NumRestricoes ; j++)
                {
                    if (i == j)
                    {
                        Base[i, j] = 1;
                    }
                    else
                    {
                        Base[i, j] = 0;
                    }
                }
            }

            // Monta vetor VarBase
            for (int i = 0 ; i < VarBase.Length ; i ++)
            {
                for (int j = NumVariaveis ; j < NumVariaveis + NumRestricoes + RestricoesErradas.Count ; j++)
                {
                    if (A[i, j] == 1)
                    {
                        VarBase[i] = j;
                    }
                }
            }

            // Monta vetor VarNBase
            VarNBase = new int[NumVariaveis + RestricoesErradas.Count];
            int k = 0;
            bool flag;
            for (int i = 0 ; i < NumVariaveis + NumRestricoes + RestricoesErradas.Count ; i++)
            {
                flag = true;
                for (int j = 0 ; j < VarBase.Length ; j++)
                {
                    if (i == VarBase[j])
                    {
                        flag = false; //Avisa que i está em VarBase
                    }
                }
                if (flag)   // se i não estava em VarBase tasca em VarNBase
                {
                    VarNBase[k] = i;
                    k++;
                }
            }

            //------------------- Fase 1 ----------------------
            if (FaseUm)
            {
                double[] FuncaoArt = new double[NumVariaveis + NumRestricoes + RestricoesErradas.Count];

                for (int i = 0; i < FuncaoArt.Length; i++)
                {
                    if (i < NumVariaveis + NumRestricoes)
                    {
                        FuncaoArt[i] = 0.0;
                    }
                    else
                    {
                        FuncaoArt[i] = 1.0;
                    }
                }

                while (true)
                {

                    //Passo 1.1 B . xB = b

                    X = new double[NumVariaveis + NumRestricoes + RestricoesErradas.Count];
                    aux = JuntarMatriz(NumRestricoes, Base, ValorRestricoes);
                    xB = gauss.CalcGauss(NumRestricoes, aux);

                    // talvez necessite arrumar por causa do retorno do Gauss
                    //Passo 1.2 X = (0, xB, 0)
                    for (int i = 0 ; i < X.Length ; i++)    // zera vetor X
                    {
                        X[i] = 0;
                    }
                    for (int i = 0 ; i < VarBase.Length ; i++) // da o valor d Xb aos que tao na base
                    {
                        X[VarBase[i]] = xB[i];
                    }

                    // Passo 2.1 bT. lambda = cB

                    bT = MatrizTransposta(NumRestricoes, Base);

                    cB = new double[VarBase.Length];

                    for (int i = 0; i < cB.Length ; i++)
                    {
                        cB[i] = FuncaoArt[VarBase[i]];
                    }

                    aux = JuntarMatriz(NumRestricoes, bT, cB);
                    Lambda = gauss.CalcGauss(NumRestricoes, aux);

                    // Passo 2.2 cN = CN - Y * aN

                    cN = new double[VarNBase.Length];

                    for (int i = 0 ; i < cN.Length ; i++)
                    {
                        double aux2 = 0;
                        for (int j = 0 ; j < Lambda.Length ; j++)
                        {
                            aux2 += Lambda[j] * A[j, VarNBase[i]];
                        }
                        cN[i] = FuncaoArt[VarNBase[i]] - aux2;
                    }

                    // Passo 2.3 Verifica cN para saber quem entra base ou se X já é solução otima

                    bool solucaoOtima = true;
                    for (int i = 0 ; i < cN.Length ; i++)
                    {
                        if (cN[i] < 0)
                        {
                            solucaoOtima = false;
                        }
                    }

                    // Se é solução otima ta errado
                    if (solucaoOtima)
                    {
                        System.Console.WriteLine("Erro: Sistema impossivel de solucionar");
                        System.Console.WriteLine("Ainda ha variaveis artificiais na base mas não ha candidatos a entrar na base.");
                        return;
                    }

                    // Pega quem entra na base e guarda em entraBase
                    double min = 0;
                    for (int i = 0 ; i < cN.Length ; i++)
                    {
                        if (cN[i] < min)
                        {
                            min = cN[i];
                            entraBase[0] = i;
                            entraBase[1] = VarNBase[i];
                        }
                    }

                    // Passo 3 By = aNk

                    nK = new double [NumRestricoes];
                    for (int i = 0 ; i < nK.Length ; i++)
                    {
                        nK[i] = A[i, entraBase[1]];
                    }

                    aux = JuntarMatriz(NumRestricoes, Base, nK);
                    y = gauss.CalcGauss(NumRestricoes, aux);

                    // Passo 4 E = min {Xb / y}  (E = SaiBase[])

                    min = 10000000; // melhorar esse min
                    flag = true;
                    for (int i = 0 ; i < y.Length ; i++)
                    {
                        if (y[i] > 0)
                        {
                            flag = false;
                            if (xB[i] / y[i] < min)
                            {
                                min = xB[i] / y[i];
                                saiBase[0] = i;
                                saiBase[1] = VarBase[i];
                            }
                        }
                    }

                    // Verifica y > 0
                    if (flag)
                    {
                        System.Console.WriteLine("Erro: Sistema sem solução");
                        System.Console.WriteLine("Não ha candidatos a sair da base porque não existe y > 0");
                        return;
                    }

                    // Passo 5 atualizar VarBase, VarNBase, Base, Verifica fim fase 1

                    VarBase[saiBase[0]] = entraBase[1];
                    VarNBase[entraBase[0]] = saiBase[1];

                    for (int i = 0 ; i < NumRestricoes ; i++)
                    {
                        Base[i, saiBase[0]] = A[i, entraBase[1]];
                    }

                    // Verifica se as variaveis Artificiais estao fora da base
                    bool fimFase1 = true;
                    for (int i = 0 ; i < VarBase.Length ; i++)
                    {
                        if (VarBase[i] >= NumVariaveis + NumRestricoes)
                        {
                            fimFase1 = false;
                        }
                    }
                    
                    // Se fimFase1
                    if (fimFase1)
                    {
                        // Remove as colunas das variaveis artificiais de A
                        for (int i = 0 ; i < RestricoesErradas.Count ; i++)
                        {
                            A = MatrizRemoveColum(NumRestricoes, tamA, A);
                            tamA--;
                        }

                        // Arruma VarNBase para a fase 2
                        VarNBase = new int[NumVariaveis];
                        k = 0;
                        for (int i = 0 ; i < NumVariaveis + NumRestricoes ; i++)
                        {
                            flag = true;
                            for (int j = 0 ; j < VarBase.Length ; j++)
                            {
                                if (i == VarBase[j])
                                {
                                    flag = false; //Avisa que i está em VarBase
                                }
                            }
                            if (flag)   // se i não estava em VarBase tasca em VarNBase
                            {
                                VarNBase[k] = i;
                                k++;
                            }
                        }

                        break;
                    }
                }

            }


            // Fase 2 --------------------------------------------------------------------------------            

            double[] FuncaoAux = new double[NumVariaveis + NumRestricoes];

            for (int i = 0; i < FuncaoAux.Length ; i++)
            {
                if (i < NumVariaveis)
                {
                    FuncaoAux[i] = Funcao[i];
                }
                else
                {
                    FuncaoAux[i] = 0.0;
                }
            }

            while (true) {
                //Passo 1.1 B . xB = b

                X = new double[NumVariaveis + NumRestricoes];
                aux = JuntarMatriz(NumRestricoes, Base, ValorRestricoes);
                xB = gauss.CalcGauss(NumRestricoes, aux);

                // talvez necessite arrumar por causa do retorno do Gauss
                //Passo 1.2 X = (0, xB, 0)
                for (int i = 0 ; i < X.Length ; i++)    // zera vetor X
                {
                    X[i] = 0;
                }

                for (int i = 0 ; i < VarBase.Length ; i++) // da o valor d Xb aos que tao na base
                {
                    X[VarBase[i]] = xB[i];
                }

                // Passo 2.1 bT. lambda = cB

                bT = MatrizTransposta(NumRestricoes, Base);

                cB = new double[VarBase.Length];

                for (int i = 0; i < cB.Length ; i++)
                {
                    cB[i] = FuncaoAux[VarBase[i]];
                }

                aux = JuntarMatriz(NumRestricoes, bT, cB);
                Lambda = gauss.CalcGauss(NumRestricoes, aux);

                // Passo 2.2 cN = CN - Y * aN

                cN = new double[VarNBase.Length];

                for (int i = 0 ; i < cN.Length ; i++)
                {
                    double aux2 = 0;
                    for (int j = 0 ; j < Lambda.Length ; j++)
                    {
                        aux2 += Lambda[j] * A[j, VarNBase[i]];
                    }
                    cN[i] = FuncaoAux[VarNBase[i]] - aux2;
                }

                // Passo 2.3 Verifica cN para saber quem entra base ou se X já é solução otima

                bool solucaoOtima = true;
                for (int i = 0 ; i < cN.Length ; i++)
                {
                    if (cN[i] < 0)
                    {
                        solucaoOtima = false;
                    }
                }

                // Se é solução otima printa resultado e break
                if (solucaoOtima)
                {
                    System.Console.WriteLine("Solução Ótima encontrada em " + it + " iterações.");
                    System.Console.Write("X = (" + X[0]);
                    for (int i = 1 ; i < X.Length ; i++)
                    {
                        System.Console.Write(" ; " + X[i]);
                    }
                    System.Console.WriteLine(")");
                    if (Min == 1)
                    {
                        System.Console.Write("Min f(x) = (");
                    }
                    else
                    {
                        System.Console.Write("Max f(x) = (");
                        for (int i = 0; i < NumVariaveis; i++)
                        {
                            Funcao[i] *= -1;
                        }
                    }
                    double resultado = 0;
                    for (int i = 0 ; i < NumVariaveis ; i++)
                    {
                        resultado += Funcao[i] * X[i];
                        System.Console.Write(Funcao[i]+ ")*(" + X[i] + ")");
                        if (i != NumVariaveis-1)
                        {
                            System.Console.Write(" + (");
                        }
                    }
                    System.Console.WriteLine();
                    if (Min == 1)
                    {
                        System.Console.WriteLine("Min f(x) = " + resultado);
                    }
                    else
                    {
                        System.Console.WriteLine("Max f(x) = " + resultado);
                    }

                    break;
                }

                // Pega quem entra na base e guarda em entraBase
                double min = 0;
                for (int i = 0 ; i < cN.Length ; i++)
                {
                    if (cN[i] < min)
                    {
                        min = cN[i];
                        entraBase[0] = i;
                        entraBase[1] = VarNBase[i];
                    }
                }

                // Passo 3 By = aNk

                nK = new double [NumRestricoes];
                for (int i = 0 ; i < nK.Length ; i++)
                {
                    nK[i] = A[i, entraBase[1]];
                }

                aux = JuntarMatriz(NumRestricoes, Base, nK);
                y = gauss.CalcGauss(NumRestricoes, aux);

                // Passo 4 E = min {Xb / y}  (E = SaiBase[])

                flag = true;
                min = 1000000000000000; // melhorar esse min
                for (int i = 0 ; i < y.Length ; i++)
                {
                    if (y[i] > 0)
                    {
                        flag = false;
                        if (xB[i] / y[i] < min)
                        {
                            min = xB[i] / y[i];
                            saiBase[0] = i;
                            saiBase[1] = VarBase[i];
                        }
                    }
                }

                // Verifica y > 0
                if (flag)
                {
                    System.Console.WriteLine("Erro: Sistema sem solução");
                    System.Console.WriteLine("Não ha candidatos a sair da base porque não existe y > 0");
                    return;
                }


                // Passo 5 atualizar VarBase, VarNBase, Base, it++

                VarBase[saiBase[0]] = entraBase[1];
                VarNBase[entraBase[0]] = saiBase[1];

                for (int i = 0 ; i < NumRestricoes ; i++)
                {
                    Base[i, saiBase[0]] = A[i, entraBase[1]];
                }

                it++;
            }
        }



        static double[,] JuntarMatriz(int tam, double[,] matriz, double[] soluc)
        {
            double[,] novaMatriz = new double[tam, tam + 1];
            for (int i = 0 ; i < tam ; i++)
            {
                for (int j = 0 ; j <= tam ; j++)
                {
                    if (j < tam)
                    {
                        novaMatriz[i, j] = matriz[i, j];
                    }
                    else
                    {
                        novaMatriz[i, j] = soluc[i];
                    }
                }
            }
            return novaMatriz;
        }

        static double[,] MatrizAddColum(int tam, int tam2,double[,] matriz, double[] soluc)
        {
            double[,] novaMatriz = new double[tam, tam2 + 1];
            for (int i = 0 ; i < tam ; i++)
            {
                for (int j = 0 ; j <= tam2 ; j++)
                {
                    if (j < tam2)
                    {
                        novaMatriz[i, j] = matriz[i, j];
                    }
                    else
                    {
                        novaMatriz[i, j] = soluc[i];
                    }
                }
            }
            return novaMatriz;
        }

        static double[,] MatrizRemoveColum(int tam, int tam2,double[,] matriz)
        {
            double[,] novaMatriz = new double[tam, tam2 - 1];
            for (int i = 0 ; i < tam ; i++)
            {
                for (int j = 0 ; j < tam2-1 ; j++)
                {
                    novaMatriz[i, j] = matriz[i, j];
                }
            }
            return novaMatriz;
        }

        static double[,] MatrizTransposta(int ordem, double[,] matrizOriginal)
        {
            double[,] matriz = new double[ordem, ordem];
            
            for (int i = 0 ; i < ordem ; i++)
            {
                for (int j = 0; j < ordem; j++)
                {
                    matriz[i, j] = matrizOriginal[j, i];
                }
            }

            return matriz;
        }


        static void printaVetor(int ordem, double[] vetor)
        {
            for (int i = 0 ; i < ordem ; i++)
            {
                System.Console.Write(vetor[i] + " ");
            }
            System.Console.WriteLine();
        }

        static void printaMatriz(int x, int y, double[,] matriz)
        {
            for (int i = 0 ; i < x ; i++)
            {
                for (int j = 0 ; j < y ; j++)
                {
                    System.Console.Write(matriz[i, j] + " ");
                }
                System.Console.WriteLine();
            }
            System.Console.WriteLine();
        }
    }
}
