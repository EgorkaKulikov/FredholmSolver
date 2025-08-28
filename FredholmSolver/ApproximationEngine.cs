using System;
using System.Collections.Generic;

namespace FredholmSolver
{
  //Алиасы типов для упрощения кода
  using Matrix = Dictionary<Tuple<int, int>, double>;
  using Vector = Dictionary<int, double>;

  public class ApproximationEngine
  {
    private Dictionary<int, double> grid;
    private Dictionary<int, double> additionalGrid;

    private Dictionary<Tuple<int, int>, double> splineProductsIntegralsCache
  = new Dictionary<Tuple<int, int>, double>();

    private SplineConstructor splineConstructor = new SplineConstructor();
    private FunctionalConstructor functionalConstructor = new FunctionalConstructor();

    public ApproximationEngine()
    {
      InitializeGrid();
      InitializeAdditionalGrid();
    }

    /// <summary>
    /// Аппроксимация в точке.
    /// </summary>
    public double ApproxInPoint(double point)
    {
      Vector coeffs = FindCoeffs(); 

      double value = 0;
      foreach (int i in coeffs.Keys)
      {
        double splineInPoint = splineConstructor.SplineValueAtPoint(i, grid, point);
        value += coeffs[i] * splineInPoint;
      }

      return value;
    }

    /// <summary>
    /// Аппроксимация в точке после применения итерации Слоана.
    /// </summary>
    public double ApproxSloanInPoint(double point)
    {
      Vector coeffs = FindCoeffs();

      double value = Configuration.F(point);
      foreach (int i in coeffs.Keys)
      {
        value += coeffs[i] * OmegaIWave(i, point);
      }

      return value;
    }
    
    public double ApproxFunction(Func<double, double> func, double point)
    {
      double value = 0;
      for (int i = -2; i <= Configuration.GridPoints - 1; i++)
      {
        value += functionalConstructor.FunctionalValue(i, func, grid, additionalGrid) *
                 splineConstructor.SplineValueAtPoint(i, grid, point) ;
      }

      return value;
    }

    /// <summary>
    /// Определение вектора коэффициентов.
    /// Решение системы уравнений (I - M)X = mu
    /// </summary>
    private Vector FindCoeffs()
    {
      Matrix matrix = IMinusMMatrix();
      Vector vector = muVector();
      Vector solution = MathUtils.SolveSystem(matrix, vector);

      return solution;
    }

    /// <summary>
    /// Построение матрицы I - M
    /// </summary>
    private Dictionary<Tuple<int, int>, double> IMinusMMatrix()
    {
      Matrix matrix = new Dictionary<Tuple<int, int>, double>();
      Matrix mMatrix = MMatrix();
      Matrix iMatrix = IMatrix();

      for (int j = -2; j <= Configuration.GridPoints - 1; j++)
      {
        for (int i = -2; i <= Configuration.GridPoints - 1; i++)
        {
          var key = Tuple.Create(j, i);
          matrix[key] = iMatrix[key] - mMatrix[key];
        }
      }

      return matrix;
    }

    /// <summary>
    /// Построение единичной матрицы I
    /// </summary>
    private Dictionary<Tuple<int, int>, double> IMatrix()
    {
      Matrix matrix = new Matrix();
      for (int j = -2; j <= Configuration.GridPoints - 1; j++)
      {
        for (int i = -2; i <= Configuration.GridPoints - 1; i++)
        {
          matrix[Tuple.Create(j, i)] = i == j ? 1 : 0;
        }
      }

      return matrix;
    }

    /// <summary>
    /// Построение матрицы M
    /// </summary>
    private Matrix MMatrix()
    {
      Matrix matrix = new Matrix();
      for (int j = -2; j <= Configuration.GridPoints - 1; j++)
      {
        for (int i = -2; i <= Configuration.GridPoints - 1; i++)
        {
          matrix[Tuple.Create(j, i)] = MMatrixElement(j, i);
        }
      }

      return matrix;
    }

    /// <summary>
    /// Построение вектора мю.
    /// </summary>
    private Dictionary<int, double> muVector()
    {
      Vector vector = new Vector();
      for (int i = -2; i <= Configuration.GridPoints - 1; i++)
      {
        vector[i] = functionalConstructor.FunctionalValue(i, Configuration.F, grid, additionalGrid);
      }

      return vector;
    }

    /// <summary>
    /// Вычисление элемента матрицы M
    /// </summary>
    private double MMatrixElement(int j, int i)
      => functionalConstructor.FunctionalValue(j, OmegaIWaveFunc(i), grid, additionalGrid);
    
    private Func<double, double> OmegaIWaveFunc(int i) => point => OmegaIWave(i, point);

    /// <summary>
    /// Вычисление Omega_i с волной
    /// Формула приведена после слов "integrals with B-spline weights"
    /// </summary>
    private double OmegaIWave(int i, double point) => OmegaWave(i)(point);
    
    private Func<double, double> OmegaWave(int i) => t =>
    {
      double FuncInT(double x) => UnderIntegralFunc(i)(t, x);
      return MathUtils.Integrate(FuncInT, Configuration.Left, Configuration.Right);
    };
    
    private Func<double, double, double> UnderIntegralFunc(int i) =>
      (t, x) => Configuration.K(t, x) * splineConstructor.SplineValueAtPoint(i, grid,x);
    

    /// <summary>
    /// Создание основной сетки.
    /// </summary>
    private void InitializeGrid()
    {
      grid = new Dictionary<int, double>();
      
      for (int i = 0; i <= Configuration.GridPoints; i++)
      {
        var theta = 1.0 * i / Configuration.GridPoints;
        
        grid[i] = Configuration.Left + theta * (Configuration.Right - Configuration.Left);
      }
      
      grid[-2] = grid[0] - 2 * Configuration.Eps;
      grid[-1] = grid[0] - Configuration.Eps;
      grid[Configuration.GridPoints + 1] = grid[Configuration.GridPoints] + Configuration.Eps;
      grid[Configuration.GridPoints + 2] = grid[Configuration.GridPoints] + 2 * Configuration.Eps;
    }
    
    /// <summary>
    /// Создание вспомогательной сетки.
    /// </summary>
    private void InitializeAdditionalGrid()
    {
      additionalGrid = new Dictionary<int, double>();
      
      for (int i = -1; i <= Configuration.GridPoints - 2; i++)
      {
        additionalGrid[i] = 0.5 * (grid[i + 1] + grid[i + 2]);
      }
      
      additionalGrid[-2] = grid[0];
      additionalGrid[Configuration.GridPoints - 1] = grid[Configuration.GridPoints]; 
    }
  }
}
