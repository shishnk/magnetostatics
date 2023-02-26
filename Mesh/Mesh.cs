namespace Magnetostatics.Mesh;

public interface IBaseMesh
{
    IReadOnlyList<Point2D> Points { get; }
    IReadOnlyList<IReadOnlyList<int>> Elements { get; }
}

public abstract class MeshBuilder
{
    protected abstract int ElementSize { get; }

    public abstract (List<Point2D>, int[][]) Build(MeshParameters meshParameters);

    protected (List<Point2D>, int[][]) BaseBuild(MeshParameters parameters)
    {
        var result = new
        {
            Points = new List<Point2D>(),
            Elements = new int[1][].Select(_ => new int[ElementSize])
                .ToArray(),
        };

        return (result.Points, result.Elements);
    }
}

public class LinearMeshBuilder : MeshBuilder
{
    protected override int ElementSize => 4;

    public override (List<Point2D>, int[][]) Build(MeshParameters parameters) => BaseBuild(parameters);
}

public class SuperMesh : IBaseMesh
{
    private readonly List<Point2D> _points;
    private readonly int[][] _elements;

    public IReadOnlyList<Point2D> Points => _points;
    public IReadOnlyList<IReadOnlyList<int>> Elements => _elements;

    public SuperMesh(MeshParameters parameters, MeshBuilder meshBuilder)
        => (_points, _elements) = meshBuilder.Build(parameters);
}

public readonly record struct Area(int Number, double Permeability, double Current, int X1, int X2, int Y1, int Y2)
{
    public static Area Parse(string line)
    {
        if (!TryParse(line, out var area))
        {
            throw new FormatException("Cant parse Area!");
        }

        return area;
    }

    public static bool TryParse(string line, out Area area)
    {
        var data = line.Split(new[] { ' ', ',' }, StringSplitOptions.RemoveEmptyEntries);

        if (data.Length != 7 || !int.TryParse(data[0], out var number) || !double.TryParse(data[1], out var mu)
            || !double.TryParse(data[2], out var current) ||
            !int.TryParse(data[3], out var x1) || !int.TryParse(data[4], out var x2) ||
            !int.TryParse(data[5], out var y1) || !int.TryParse(data[6], out var y2))
        {
            area = default;
            return false;
        }

        area = new(number, mu, current, x1, x2, y1, y2);
        return true;
    }
}

public class MeshParameters
{
    public required ImmutableArray<double> LinesX { get; init; }
    public required ImmutableArray<double> LinesY { get; init; }
    public required ImmutableArray<int> SplitsX { get; init; }
    public required ImmutableArray<int> SplitsY { get; init; }
    public required ImmutableArray<int> Kx { get; init; }
    public required ImmutableArray<int> Ky { get; init; }
    public required (int, int) Nesting { get; init; }
    public required ImmutableArray<Area> Areas { get; init; }

    public MeshParameters(IEnumerable<double> linesX, IEnumerable<double> linesY, IEnumerable<int> splitsX,
        IEnumerable<int> splitsY, IEnumerable<int> kx, IEnumerable<int> ky, (int, int) nesting,
        IEnumerable<Area> areas)
    {
        LinesX = linesX.ToImmutableArray();
        LinesY = linesY.ToImmutableArray();
        SplitsX = splitsX.ToImmutableArray();
        SplitsY = splitsY.ToImmutableArray();
        Kx = kx.ToImmutableArray();
        Ky = ky.ToImmutableArray();
        Nesting = nesting;
        Areas = areas.ToImmutableArray();
    }


    public MeshParameters(string path)
    {
        if (!File.Exists(path))
        {
            throw new("File does not exist");
        }

        using var sr = new StreamReader(path);
        LinesX = sr.ReadLine()!.Split().Select(double.Parse).ToImmutableArray();
        LinesY = sr.ReadLine()!.Split().Select(double.Parse).ToImmutableArray();
        SplitsX = sr.ReadLine()!.Split().Select(int.Parse).ToImmutableArray();
        SplitsY = sr.ReadLine()!.Split().Select(int.Parse).ToImmutableArray();
        Nesting = sr.ReadLine()!.Split().Select(line => line.Split())
            .Select(value => (int.Parse(value[0]), int.Parse(value[1]))).DefaultIfEmpty((0, 0)).FirstOrDefault();
        Areas = sr.ReadToEnd().Split("\n").Select(Area.Parse).ToImmutableArray();
    }
}