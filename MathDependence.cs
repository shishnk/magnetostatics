namespace Magnetostatics;

public class MathDependence
{
    public string Argument { get; }
    public string Function { get; }
    public (double Function, double Argument)[]? Data { get; private set; }

    public MathDependence(string function, string argument)
        => (Function, Argument) = (function, argument);

    public void LoadData(string path)
    {
        if (!File.Exists(path)) throw new("File does not exist");

        using var sr = new StreamReader(path);

        Data = sr.ReadToEnd().Split("\n").Select(line => line.Split(" ", StringSplitOptions.RemoveEmptyEntries))
            .Select(line => (double.Parse(line[0]), double.Parse(line[1]))).ToArray();
    }

    public override string ToString() => $"Dependence is {Function}({Argument})";
}