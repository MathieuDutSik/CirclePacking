TheDodec:=ArchimedeanPolyhedra("Dodecahedron");

PLwork:=GoldbergConstruction(TheDodec, 2, 1);

TestFile:="TheRes21.out";
OutputToMoharProgram(TestFile, PLwork);

