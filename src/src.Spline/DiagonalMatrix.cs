namespace Magnetostatics.src.Spline;

/// <summary>
/// Class implement seven diagonal matrix without zero diagonals
/// </summary>
public class SevenDiagonalMatrix
{
    public int Size { get; }
    public double[][] Diags { get; }

    public double this[int i, int j]
    {
        get => Diags[i][j];
        set => Diags[i][j] = value;
    }

    public SevenDiagonalMatrix(int mainDiagonalSize)
    {
        Size = mainDiagonalSize;
        Diags = new double[7][];
        Diags[0] = new double[mainDiagonalSize];
        Diags[1] = new double[mainDiagonalSize - 1];
        Diags[4] = new double[Diags[1].Length];
        Diags[2] = new double[mainDiagonalSize - 2];
        Diags[5] = new double[Diags[2].Length];
        Diags[3] = new double[mainDiagonalSize - 3];
        Diags[6] = new double[Diags[3].Length];
    }
}