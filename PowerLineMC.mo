package PowerLineMC
  model CompareHypTStep "Compares different line models using sine wave source"
    parameter Real Rload = 1000;
    //load resistance value
    parameter Real r1 = 1e-6, g1 = 1e-12, l1 = 1e-6, c1 = 1e-11, len = 100e3;
    //resistance, conductance, inductance, capacitance per unit length, length
    extends Modelica.Icons.Example;
    import Modelica.ComplexMath.j;
    import Modelica.ComplexMath.real;
    import Modelica.ComplexMath.imag;
    import Modelica.ComplexMath.sinh;
    import Modelica.ComplexMath.tanh;
    import Modelica.Constants.pi;
    final parameter Real w = 2*pi*50;
    parameter Real c = 1/sqrt(l1*c1), z0 = sqrt(l1/c1), td = len/c;
    parameter Complex Z0compl = Modelica.ComplexMath.sqrt((r1 + j*w*l1)/(j*w*c1));
    parameter Complex K = Modelica.ComplexMath.sqrt((r1 + j*w*l1)*j*w*c1);
    parameter Real rHypLong = real(Z0compl*tanh(K*len/2));
    parameter Real lHypLong = imag(Z0compl*tanh(K*len/2))/w;
    parameter Real rHypTransv = real(sinh(K*len)/Z0compl);
    parameter Real cHypTransv = imag(sinh(K*len)/Z0compl)/w;
    Modelica.Electrical.Analog.Lines.TLine1 tLine(Z0 = z0, TD = td) annotation(
      Placement(visible = true, transformation(origin = {-2, -16}, extent = {{-30, -46}, {-10, -26}}, rotation = 0)));
    Modelica.Electrical.Analog.Sources.SignalVoltage sourceV annotation(
      Placement(visible = true, transformation(origin = {-54, -38}, extent = {{-10, 10}, {10, -10}}, rotation = -90)));
    Modelica.Electrical.Analog.Basic.Ground ground annotation(
      Placement(visible = true, transformation(origin = {-2, -6}, extent = {{-62, -76}, {-42, -56}}, rotation = 0)));
    Modelica.Electrical.Analog.Basic.Resistor Rdistr(R = Rload) annotation(
      Placement(visible = true, transformation(origin = {2, -46}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Electrical.Analog.Lines.OLine oLine10(N = 10, r = r1, l = l1, g = g1, c = c1, length = len) annotation(
      Placement(visible = true, transformation(origin = {58, 28}, extent = {{26, 20}, {46, 40}}, rotation = 0)));
    Modelica.Electrical.Analog.Basic.Resistor R10(R = Rload) annotation(
      Placement(visible = true, transformation(origin = {122, 48}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Electrical.Analog.Basic.Ground gTr annotation(
      Placement(visible = true, transformation(origin = {-2, -6}, extent = {{-14, -76}, {6, -56}}, rotation = 0)));
    Modelica.Electrical.Analog.Lines.OLine oLine50(r = r1, l = l1, g = g1, c = c1, length = len, N = 50) annotation(
      Placement(visible = true, transformation(origin = {20, -18}, extent = {{26, -42}, {46, -22}}, rotation = 0)));
    Modelica.Electrical.Analog.Basic.Resistor R50(R = Rload) annotation(
      Placement(visible = true, transformation(origin = {86, -52}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Electrical.Analog.Basic.Capacitor CHypTransv(C = cHypTransv) annotation(
      Placement(visible = true, transformation(origin = {8, 24}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Electrical.Analog.Basic.Resistor RHypLoad(R = Rload) annotation(
      Placement(visible = true, transformation(origin = {75, 23}, extent = {{-11, -11}, {11, 11}}, rotation = -90)));
    Modelica.Electrical.Analog.Basic.Resistor RHypLong(R = rHypLong) annotation(
      Placement(visible = true, transformation(origin = {-32, 44}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Electrical.Analog.Basic.Inductor LHypLong(L = lHypLong) annotation(
      Placement(visible = true, transformation(origin = {-8, 44}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Electrical.Analog.Basic.Ground gHyp annotation(
      Placement(visible = true, transformation(origin = {96, 44}, extent = {{-14, -76}, {6, -56}}, rotation = 0)));
    Modelica.Electrical.Analog.Basic.Resistor RHypTransv(R = rHypTransv) annotation(
      Placement(visible = true, transformation(origin = {9, -3}, extent = {{-11, -11}, {11, 11}}, rotation = -90)));
    Modelica.Blocks.Sources.Step step(startTime = 200e-6) annotation(
      Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = -90, origin = {-80, 30})));
    Modelica.Blocks.Continuous.FirstOrder firstOrder(T = 20e-6) annotation(
      Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = -90, origin = {-80, -6})));
    Modelica.Electrical.Analog.Basic.Inductor lHypLong1(L = lHypLong) annotation(
      Placement(visible = true, transformation(origin = {54, 44}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Electrical.Analog.Basic.Resistor rHypLong1(R = rHypLong) annotation(
      Placement(visible = true, transformation(origin = {26, 44}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(ground.p, sourceV.n) annotation(
      Line(points = {{-54, -62}, {-54, -48}}, color = {0, 0, 255}));
    connect(Rdistr.p, tLine.p2) annotation(
      Line(points = {{2, -36}, {-12, -36}, {-12, -42}}, color = {0, 0, 255}));
    connect(Rdistr.n, tLine.n2) annotation(
      Line(points = {{2, -56}, {2, -62}, {-12, -62}}, color = {0, 0, 255}));
    connect(R10.p, oLine10.p2) annotation(
      Line(points = {{122, 58}, {104, 58}}, color = {0, 0, 255}));
    connect(R10.n, oLine10.p3) annotation(
      Line(points = {{122, 38}, {94, 38}, {94, 48}}, color = {0, 0, 255}));
    connect(gTr.p, tLine.n2) annotation(
      Line(points = {{-6, -62}, {-12, -62}}, color = {0, 0, 255}));
    connect(R50.p, oLine50.p2) annotation(
      Line(points = {{86, -42}, {76, -42}, {76, -50}, {66, -50}}, color = {0, 0, 255}));
    connect(R50.n, oLine50.p3) annotation(
      Line(points = {{86, -62}, {56, -62}, {56, -60}}, color = {0, 0, 255}));
    connect(oLine10.p1, sourceV.p) annotation(
      Line(points = {{84, 58}, {-54, 58}, {-54, -28}}, color = {0, 0, 255}));
    connect(oLine50.p1, sourceV.p) annotation(
      Line(points = {{46, -50}, {20, -50}, {20, -30}, {-32, -30}, {-32, -26}, {-54, -26}, {-54, -28}}, color = {0, 0, 255}));
    connect(tLine.p1, sourceV.p) annotation(
      Line(points = {{-32, -42}, {-32, -26}, {-54, -26}, {-54, -28}}, color = {0, 0, 255}));
    connect(tLine.n1, ground.p) annotation(
      Line(points = {{-32, -62}, {-54, -62}}, color = {0, 0, 255}));
    connect(RHypLong.n, LHypLong.p) annotation(
      Line(points = {{-22, 44}, {-18, 44}}, color = {0, 0, 255}));
    connect(RHypLong.p, sourceV.p) annotation(
      Line(points = {{-42, 44}, {-54, 44}, {-54, -28}}, color = {0, 0, 255}));
    connect(gHyp.p, RHypLoad.n) annotation(
      Line(points = {{92, -12}, {92, -2}, {75, -2}, {75, 12}}, color = {0, 0, 255}));
    connect(RHypTransv.p, CHypTransv.n) annotation(
      Line(points = {{9, 8}, {8, 8}, {8, 14}}, color = {0, 0, 255}));
    connect(R10.n, gHyp.p) annotation(
      Line(points = {{122, 38}, {122, -2}, {92, -2}, {92, -12}}, color = {0, 0, 255}));
    connect(R50.n, RHypLoad.n) annotation(
      Line(points = {{86, -62}, {106, -62}, {106, -2}, {75, -2}, {75, 12}}, color = {0, 0, 255}));
    connect(firstOrder.u, step.y) annotation(
      Line(points = {{-80, 6}, {-80, 19}}, color = {0, 0, 127}, smooth = Smooth.None));
    connect(firstOrder.y, sourceV.v) annotation(
      Line(points = {{-80, -17}, {-80, -38}, {-66, -38}}, color = {0, 0, 127}));
    connect(LHypLong.n, rHypLong1.p) annotation(
      Line(points = {{2, 44}, {16, 44}}, color = {0, 0, 255}));
    connect(lHypLong1.n, RHypLoad.p) annotation(
      Line(points = {{64, 44}, {75, 44}, {75, 34}}, color = {0, 0, 255}));
    connect(CHypTransv.p, LHypLong.n) annotation(
      Line(points = {{8, 34}, {8, 44}, {2, 44}}, color = {0, 0, 255}));
    connect(RHypTransv.n, RHypLoad.n) annotation(
      Line(points = {{9, -14}, {8, -14}, {8, -18}, {75, -18}, {75, 12}}, color = {0, 0, 255}));
    connect(rHypLong1.n, lHypLong1.p) annotation(
      Line(points = {{36, 44}, {44, 44}}, color = {0, 0, 255}));
    annotation(
      Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, 80}, {140, -80}})),
      Icon(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = false)),
      experiment(StopTime = 0.003),
      __Dymola_experimentSetupOutput,
      Documentation(info = "<html><head></head><body>Aggiunge a smoothStep il confronto col pi-greco con formule iperboliche.<div>Per ora però la sollecitazione è solo sinusoidale</div></body></html>"));
  end CompareHypTStep;

  model CompareHypTSine "Compares different line models using sine wave source"
    parameter Real Rload = 1000;
    //load resistance value
    parameter Real r1 = 1e-6, g1 = 1e-12, l1 = 1e-6, c1 = 1e-11, len = 500e3;
    //resistance, conductance, inductance, capacitance per unit length, length
    extends Modelica.Icons.Example;
    import Modelica.ComplexMath.j;
    import Modelica.ComplexMath.real;
    import Modelica.ComplexMath.imag;
    import Modelica.ComplexMath.sinh;
    import Modelica.ComplexMath.tanh;
    import Modelica.Constants.pi;
    final parameter Real w = 2*pi*50;
    parameter Real c = 1/sqrt(l1*c1), z0 = sqrt(l1/c1), td = len/c;
    parameter Complex Z0compl = Modelica.ComplexMath.sqrt((r1 + j*w*l1)/(j*w*c1));
    parameter Complex K = Modelica.ComplexMath.sqrt((r1 + j*w*l1)*j*w*c1);
    parameter Real rHypLong = real(Z0compl*tanh(K*len/2));
    parameter Real lHypLong = imag(Z0compl*tanh(K*len/2))/w;
    parameter Real rHypTransv = real(sinh(K*len)/Z0compl);
    parameter Real cHypTransv = imag(sinh(K*len)/Z0compl)/w;
    Modelica.Electrical.Analog.Lines.TLine1 tLine(Z0 = z0, TD = td) annotation(
      Placement(visible = true, transformation(origin = {-2, -16}, extent = {{-30, -46}, {-10, -26}}, rotation = 0)));
    Modelica.Electrical.Analog.Sources.SignalVoltage sourceV annotation(
      Placement(visible = true, transformation(origin = {-54, -38}, extent = {{-10, 10}, {10, -10}}, rotation = -90)));
    Modelica.Electrical.Analog.Basic.Ground ground annotation(
      Placement(visible = true, transformation(origin = {-2, -6}, extent = {{-62, -76}, {-42, -56}}, rotation = 0)));
    Modelica.Electrical.Analog.Basic.Resistor Rdistr(R = Rload) annotation(
      Placement(visible = true, transformation(origin = {2, -46}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Electrical.Analog.Lines.OLine oLine1(N = 1, c = c1, g = g1, l = l1, length = len, r = r1) annotation(
      Placement(visible = true, transformation(origin = {58, 28}, extent = {{26, 20}, {46, 40}}, rotation = 0)));
    Modelica.Electrical.Analog.Basic.Resistor R1(R = Rload) annotation(
      Placement(visible = true, transformation(origin = {122, 48}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Electrical.Analog.Basic.Ground gTr annotation(
      Placement(visible = true, transformation(origin = {-2, -6}, extent = {{-14, -76}, {6, -56}}, rotation = 0)));
    Modelica.Electrical.Analog.Lines.OLine oLine10(r = r1, l = l1, g = g1, c = c1, length = len, N = 10) annotation(
      Placement(visible = true, transformation(origin = {20, -18}, extent = {{26, -42}, {46, -22}}, rotation = 0)));
    Modelica.Electrical.Analog.Basic.Resistor R10(R = Rload) annotation(
      Placement(visible = true, transformation(origin = {86, -52}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Electrical.Analog.Basic.Capacitor CHypTransv(C = cHypTransv) annotation(
      Placement(visible = true, transformation(origin = {8, 24}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Electrical.Analog.Basic.Resistor RHypLoad(R = Rload) annotation(
      Placement(visible = true, transformation(origin = {75, 23}, extent = {{-11, -11}, {11, 11}}, rotation = -90)));
    Modelica.Electrical.Analog.Basic.Resistor RHypLong(R = rHypLong) annotation(
      Placement(visible = true, transformation(origin = {-32, 44}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Electrical.Analog.Basic.Inductor LHypLong(L = lHypLong) annotation(
      Placement(visible = true, transformation(origin = {-8, 44}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Electrical.Analog.Basic.Ground gHyp annotation(
      Placement(visible = true, transformation(origin = {96, 44}, extent = {{-14, -76}, {6, -56}}, rotation = 0)));
    Modelica.Electrical.Analog.Basic.Resistor RHypTransv(R = rHypTransv) annotation(
      Placement(visible = true, transformation(origin = {9, -3}, extent = {{-11, -11}, {11, 11}}, rotation = -90)));
    Modelica.Blocks.Continuous.FirstOrder firstOrder(T = 20e-6) annotation(
      Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = -90, origin = {-80, -6})));
    Modelica.Electrical.Analog.Basic.Inductor lHypLong1(L = lHypLong) annotation(
      Placement(visible = true, transformation(origin = {54, 44}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Electrical.Analog.Basic.Resistor rHypLong1(R = rHypLong) annotation(
      Placement(visible = true, transformation(origin = {26, 44}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Sine sine(f = 50) annotation(
      Placement(visible = true, transformation(origin = {-80, 30}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  equation
    connect(ground.p, sourceV.n) annotation(
      Line(points = {{-54, -62}, {-54, -48}}, color = {0, 0, 255}));
    connect(Rdistr.p, tLine.p2) annotation(
      Line(points = {{2, -36}, {-12, -36}, {-12, -42}}, color = {0, 0, 255}));
    connect(Rdistr.n, tLine.n2) annotation(
      Line(points = {{2, -56}, {2, -62}, {-12, -62}}, color = {0, 0, 255}));
    connect(R1.p, oLine1.p2) annotation(
      Line(points = {{122, 58}, {104, 58}}, color = {0, 0, 255}));
    connect(R1.n, oLine1.p3) annotation(
      Line(points = {{122, 38}, {94, 38}, {94, 48}}, color = {0, 0, 255}));
    connect(gTr.p, tLine.n2) annotation(
      Line(points = {{-6, -62}, {-12, -62}}, color = {0, 0, 255}));
    connect(R10.p, oLine10.p2) annotation(
      Line(points = {{86, -42}, {76, -42}, {76, -50}, {66, -50}}, color = {0, 0, 255}));
    connect(R10.n, oLine10.p3) annotation(
      Line(points = {{86, -62}, {56, -62}, {56, -60}}, color = {0, 0, 255}));
    connect(oLine1.p1, sourceV.p) annotation(
      Line(points = {{84, 58}, {-54, 58}, {-54, -28}}, color = {0, 0, 255}));
    connect(oLine10.p1, sourceV.p) annotation(
      Line(points = {{46, -50}, {20, -50}, {20, -30}, {-32, -30}, {-32, -26}, {-54, -26}, {-54, -28}}, color = {0, 0, 255}));
    connect(tLine.p1, sourceV.p) annotation(
      Line(points = {{-32, -42}, {-32, -26}, {-54, -26}, {-54, -28}}, color = {0, 0, 255}));
    connect(tLine.n1, ground.p) annotation(
      Line(points = {{-32, -62}, {-54, -62}}, color = {0, 0, 255}));
    connect(RHypLong.n, LHypLong.p) annotation(
      Line(points = {{-22, 44}, {-18, 44}}, color = {0, 0, 255}));
    connect(RHypLong.p, sourceV.p) annotation(
      Line(points = {{-42, 44}, {-54, 44}, {-54, -28}}, color = {0, 0, 255}));
    connect(gHyp.p, RHypLoad.n) annotation(
      Line(points = {{92, -12}, {92, -2}, {75, -2}, {75, 12}}, color = {0, 0, 255}));
    connect(RHypTransv.p, CHypTransv.n) annotation(
      Line(points = {{9, 8}, {8, 8}, {8, 14}}, color = {0, 0, 255}));
    connect(R1.n, gHyp.p) annotation(
      Line(points = {{122, 38}, {122, -2}, {92, -2}, {92, -12}}, color = {0, 0, 255}));
    connect(R10.n, RHypLoad.n) annotation(
      Line(points = {{86, -62}, {106, -62}, {106, -2}, {75, -2}, {75, 12}}, color = {0, 0, 255}));
    connect(firstOrder.y, sourceV.v) annotation(
      Line(points = {{-80, -17}, {-80, -38}, {-66, -38}}, color = {0, 0, 127}));
    connect(LHypLong.n, rHypLong1.p) annotation(
      Line(points = {{2, 44}, {16, 44}}, color = {0, 0, 255}));
    connect(lHypLong1.n, RHypLoad.p) annotation(
      Line(points = {{64, 44}, {75, 44}, {75, 34}}, color = {0, 0, 255}));
    connect(CHypTransv.p, LHypLong.n) annotation(
      Line(points = {{8, 34}, {8, 44}, {2, 44}}, color = {0, 0, 255}));
    connect(RHypTransv.n, RHypLoad.n) annotation(
      Line(points = {{9, -14}, {8, -14}, {8, -18}, {75, -18}, {75, 12}}, color = {0, 0, 255}));
    connect(rHypLong1.n, lHypLong1.p) annotation(
      Line(points = {{36, 44}, {44, 44}}, color = {0, 0, 255}));
    connect(firstOrder.u, sine.y) annotation(
      Line(points = {{-80, 6}, {-80, 20}}, color = {0, 0, 127}));
    annotation(
      Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-100, 80}, {140, -80}})),
      Icon(coordinateSystem(extent = {{-100, 100}, {100, -100}}, preserveAspectRatio = false)),
      experiment(StopTime = 0.09, StartTime = 0, Tolerance = 1e-06, Interval = 4.5e-05),
      __Dymola_experimentSetupOutput,
      Documentation(info = "<html><head></head><body>Aggiunge a smoothStep il confronto col pi-greco con formule iperboliche.<div>Per ora però la sollecitazione è solo sinusoidale</div></body></html>"));
  end CompareHypTSine;

  model PowerLineWithFence
    extends Modelica.Icons.Example;
    parameter Real Rl = 10000, Rb = 1e6;
    Modelica.Electrical.Analog.Lines.M_OLine line(N = 2, c = Ccomp, g = fill(1e-12, div(g.n*(g.n + 1), 2)), l = Lcomp, length(displayUnit = "km") = 100000, lines = g.n, r = g.R1) annotation(
      Placement(visible = true, transformation(origin = {-16, 8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Electrical.Analog.Sources.SineVoltage v2(V = v1.V, phase = -2.0943951023932, f = v1.f) annotation(
      Placement(visible = true, transformation(origin = {-56, 8}, extent = {{10, 10}, {-10, -10}}, rotation = 180)));
    Modelica.Electrical.Analog.Basic.Resistor rgUL(R = 1) annotation(
      Placement(visible = true, transformation(origin = {-76, -32}, extent = {{-10, 10}, {10, -10}}, rotation = 270)));
    Modelica.Electrical.Analog.Basic.Ground ground annotation(
      Placement(visible = true, transformation(origin = {-6, -2}, extent = {{-80, -70}, {-60, -50}}, rotation = 0)));
    Modelica.Electrical.Analog.Sources.SineVoltage v3(V = v1.V, phase = 2.0943951023932, f = v1.f) annotation(
      Placement(visible = true, transformation(origin = {-56, -12}, extent = {{10, 10}, {-10, -10}}, rotation = 180)));
    Modelica.Electrical.Analog.Basic.Resistor r1(R = Rl) annotation(
      Placement(visible = true, transformation(origin = {24, 28}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Electrical.Analog.Basic.Resistor r2(R = r1.R) annotation(
      Placement(visible = true, transformation(origin = {24, 8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Electrical.Analog.Basic.Resistor r3(R = r1.R) annotation(
      Placement(visible = true, transformation(origin = {24, -12}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  parameter PowerLineMC.Functions.LineGeometry g(n = 4, x = {0, -3.048, 3.048, -9.144}, y = {12.192, 12.192, 12.192, 3.048}, r = 1e-3/2*{12.7, 12.7, 12.7, 4.064}, R1 = 1e-3*{0.348, 0.348, 0.348, 1.802}, k_s = {0.287, 0.287, 0.287, 0.779}, f = 60) annotation(
      Placement(visible = true, transformation(origin = {-76, -134}, extent = {{40, 60}, {60, 80}}, rotation = 0)));
    Modelica.Electrical.Analog.Sources.SineVoltage v1(V = 345e3*sqrt(2/3), f = 60) annotation(
      Placement(visible = true, transformation(origin = {-56, 28}, extent = {{10, 10}, {-10, -10}}, rotation = 180)));
    Modelica.Electrical.Analog.Basic.Resistor body(R = Rb) annotation(
      Placement(visible = true, transformation(origin = {-6, -32}, extent = {{-10, 10}, {10, -10}}, rotation = 270)));
    parameter Real Rcomp[div(g.n*(g.n + 1), 2)](each final unit = "Ohm/m", each fixed = false) "Compact resistance matrix";
    parameter Real Xcomp[div(g.n*(g.n + 1), 2)](each final unit = "Ohm/m", each fixed = false) "Compact reactance matrix";
    parameter Real Lcomp[div(g.n*(g.n + 1), 2)](each final unit = "H/m", each fixed = false) "Compact inductance";
    parameter Real Ccomp[div(g.n*(g.n + 1), 2)](each final unit = "F/m") = Functions.lineCmatrix(n = g.n, x = g.x, y = g.y, r = g.r);
  initial algorithm
    (Rcomp, Xcomp, Lcomp) := Functions.lineZmatrix(n = g.n, x = g.x, y = g.y, r = g.r, R1 = g.R1, k_s = g.k_s, rho = g.rho, f = g.f);
  equation
    connect(ground.p, rgUL.n) annotation(
      Line(points = {{-76, -52}, {-76, -42}}, color = {0, 0, 255}, pattern = LinePattern.Solid));
    connect(r1.p, line.n[1]) annotation(
      Line(points = {{14, 28}, {-6, 28}, {-6, 8}}, color = {0, 0, 255}, pattern = LinePattern.Solid));
    connect(r2.p, line.n[2]) annotation(
      Line(points = {{14, 8}, {-6, 8}}, color = {0, 0, 255}, pattern = LinePattern.Solid));
    connect(r1.n, r2.n) annotation(
      Line(points = {{34, 28}, {44, 28}, {44, 8}, {34, 8}}, color = {0, 0, 255}, pattern = LinePattern.Solid));
    connect(r2.n, r3.n) annotation(
      Line(points = {{34, 8}, {44, 8}, {44, -12}, {34, -12}}, color = {0, 0, 255}, pattern = LinePattern.Solid));
    connect(r3.p, line.n[3]) annotation(
      Line(points = {{14, -12}, {-6, -12}, {-6, 8}}, color = {0, 0, 255}, pattern = LinePattern.Solid));
    connect(v1.n, line.p[1]) annotation(
      Line(points = {{-46, 28}, {-26, 28}, {-26, 8}}, color = {0, 0, 255}, pattern = LinePattern.Solid));
    connect(v2.n, line.p[2]) annotation(
      Line(points = {{-46, 8}, {-26, 8}}, color = {0, 0, 255}, pattern = LinePattern.Solid));
    connect(v3.n, line.p[3]) annotation(
      Line(points = {{-46, -12}, {-26, -12}, {-26, 8}}, color = {0, 0, 255}, pattern = LinePattern.Solid));
    connect(v1.p, rgUL.p) annotation(
      Line(points = {{-66, 28}, {-76, 28}, {-76, -22}}, color = {0, 0, 255}, pattern = LinePattern.Solid));
    connect(v2.p, rgUL.p) annotation(
      Line(points = {{-66, 8}, {-76, 8}, {-76, -22}}, color = {0, 0, 255}, pattern = LinePattern.Solid));
    connect(v3.p, rgUL.p) annotation(
      Line(points = {{-66, -12}, {-76, -12}, {-76, -22}}, color = {0, 0, 255}, pattern = LinePattern.Solid));
    connect(body.p, line.n[4]) annotation(
      Line(points = {{-6, -22}, {-6, 8}}, color = {0, 0, 255}, pattern = LinePattern.Solid));
    connect(body.n, ground.p) annotation(
      Line(points = {{-6, -42}, {-6, -52}, {-76, -52}}, color = {0, 0, 255}, pattern = LinePattern.Solid));
    connect(body.n, r3.n) annotation(
      Line(points = {{-6, -42}, {-6, -52}, {44, -52}, {44, -12}, {34, -12}}, color = {0, 0, 255}, pattern = LinePattern.Solid));
    annotation(
      Diagram(coordinateSystem(extent = {{-100, -80}, {60, 40}})),
      Documentation(info = "<html>
<p>
This example shows the usage of <a href=\"modelica://Modelica.Electrical.Analog.Lines.M_OLine\">M_OLine</a>
in a overhead PowerLine. The line geometry is the one shown in figure 4.11 of
[<a href=\"modelica://Modelica.Electrical.Analog.UsersGuide.References\">Theory Book</a>].
</p>
<p>
A common rule-of-thumb states that for steady-state analysis at industrial frequency each
segment (PI model) should be not longer than one tenth of the wave length. For 50&nbsp;Hz,
the wavelength is 6000&nbsp;km, therefore each pi should not overcome 60&nbsp;km. This is
why for our 100&nbsp;km line we have chosen just two trunks.
</p>
<p>
This example has some limitations:
</p>
<ul>
  <li>
    the results are valid only for their steady-state trend. If we want to simulate transients
    effectively, the number of trunks must be enlarged consistently. See example
    &quot;CompareLineTrunks&quot;.
  </li>
  <li>
    since it is based on <a href=\"modelica://Modelica.Electrical.Analog.Lines.M_OLine\">M_OLine</a>,
    which requires zero-valued resistances on the off-diagonal elements of Z impedance matrix,
    it does not take into account the influence of the soil resistance into Z matrix; however,
    for simulations like this one in which the current through the soil is negligible and we
    analyse the results at power frequency where soil resistance is very low, this implies
    negligible errors.
  </li>
</ul>
<p>
The considered line is a three-phase line with a fence all along the line (conductor #4).
An interesting result of this example is the fence steady-state voltage due to the capacitive
coupling with the other wires.
</p>

<h5>Test 1:</h5>
<p>
Simulate the model as it is: a person (body) is supposed to be touching the fence, and
absorbs the current <code>body.i</code>, around 1&nbsp;A (RMS), by large lethal.
</p>

<h5>Test 2:</h5>
<p>
Simulate the model changing the body resistance to 1e6&nbsp;&ohm; or removing body: the steady-state
fence voltage <code>line.n[4].v</code> has changed from 1429&nbsp;V to 5594&nbsp;V (peak).
So there exists a marked difference between the open-circuit voltage and actual body
voltage which, although being much lower, is still very dangerous.
</p>

<h5>Test 3:</h5>
<p>
Simulate the original model reducing the line length to 3&nbsp;km. The body current has fallen
to around 30&nbsp;mA (RMS), which is considered the maximum safe current a body can resist over
5&nbsp;s.
</p>
<p>
Note that the fence voltage before the human contact is roughtly the same for all lengths while
the current flowing when contact occurs is very dependent on the line length.
</p>
</html>", revisions = "<html>
<p>June, 2023 Massimo Ceraolo of the University of Pisa </p>
<p>originally created </p>
</html>"),
      experiment(StopTime = 0.04, Interval = 1e-05, Tolerance = 1e-06),
      Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}})));
  end PowerLineWithFence;

  model PowerLineWithFenceSC
    extends Modelica.Icons.Example;
    Modelica.Electrical.Analog.Lines.M_OLine line(N = 2, c = Ccomp, g = fill(1e-12, div(g.n*(g.n + 1), 2)), l = Lcomp, length(displayUnit = "km") = 2000, lines = g.n, r = g.R1) annotation(
      Placement(visible = true, transformation(origin = {-10, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Electrical.Analog.Sources.SineVoltage v2(V = v1.V, phase = -2.0943951023932, f = v1.f) annotation(
      Placement(visible = true, transformation(origin = {-50, 10}, extent = {{10, 10}, {-10, -10}}, rotation = 180)));
    Modelica.Electrical.Analog.Basic.Resistor rgUL(R = 1) annotation(
      Placement(visible = true, transformation(origin = {-70, -30}, extent = {{-10, 10}, {10, -10}}, rotation = 270)));
    Modelica.Electrical.Analog.Basic.Ground ground annotation(
      Placement(visible = true, transformation(origin = {0, 0}, extent = {{-80, -70}, {-60, -50}}, rotation = 0)));
    Modelica.Electrical.Analog.Sources.SineVoltage v3(V = v1.V, phase = 2.0943951023932, f = v1.f) annotation(
      Placement(visible = true, transformation(origin = {-50, -10}, extent = {{10, 10}, {-10, -10}}, rotation = 180)));
    Modelica.Electrical.Analog.Basic.Resistor r1(R = 1000) annotation(
      Placement(visible = true, transformation(origin = {30, 30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Electrical.Analog.Basic.Resistor r2(R = 1e9) annotation(
      Placement(visible = true, transformation(origin = {30, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Electrical.Analog.Basic.Resistor r3(R = 1e9) annotation(
      Placement(visible = true, transformation(origin = {30, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Functions.LineGeometry g(n = 4, x = {0, -3.048, 3.048, -9.144}, y = {12.192, 12.192, 12.192, 3.048}, r = 1e-3/2*{12.7, 12.7, 12.7, 4.064}, R1 = 1e-3*{0.348, 0.348, 0.348, 1.802}, k_s = {0.287, 0.287, 0.287, 0.779}, f = 60) annotation(
      Placement(visible = true, transformation(origin = {-56, -132}, extent = {{40, 60}, {60, 80}}, rotation = 0)));
    Modelica.Electrical.Analog.Sources.SineVoltage v1(V = 345e3*sqrt(2/3), f = 60) annotation(
      Placement(visible = true, transformation(origin = {-50, 30}, extent = {{10, 10}, {-10, -10}}, rotation = 180)));
    Modelica.Electrical.Analog.Basic.Resistor body(R = 1e-3) annotation(
      Placement(visible = true, transformation(origin = {0, -30}, extent = {{-10, 10}, {10, -10}}, rotation = 270)));
    parameter Real Rcomp[div(g.n*(g.n + 1), 2)](each final unit = "Ohm/m", each fixed = false) "Compact resistance matrix";
    parameter Real Xcomp[div(g.n*(g.n + 1), 2)](each final unit = "Ohm/m", each fixed = false) "Compact reactance matrix";
    parameter Real Lcomp[div(g.n*(g.n + 1), 2)](each final unit = "H/m", each fixed = false) "Compact inductance";
    parameter Real Ccomp[div(g.n*(g.n + 1), 2)](each final unit = "F/m") = Functions.lineCmatrix(n = g.n, x = g.x, y = g.y, r = g.r);
    Modelica.Electrical.Analog.Basic.Resistor Conn1(R = 1e-3) annotation(
      Placement(visible = true, transformation(origin = {-20, -30}, extent = {{-10, 10}, {10, -10}}, rotation = 270)));
  initial algorithm
    (Rcomp, Xcomp, Lcomp) := Functions.lineZmatrix(n = g.n, x = g.x, y = g.y, r = g.r, R1 = g.R1, k_s = g.k_s, rho = g.rho, f = g.f);
  equation
    connect(ground.p, rgUL.n) annotation(
      Line(points = {{-70, -50}, {-70, -40}}, color = {0, 0, 255}, pattern = LinePattern.Solid));
    connect(r1.p, line.n[1]) annotation(
      Line(points = {{20, 30}, {0, 30}, {0, 10}}, color = {0, 0, 255}, pattern = LinePattern.Solid));
    connect(r2.p, line.n[2]) annotation(
      Line(points = {{20, 10}, {0, 10}}, color = {0, 0, 255}, pattern = LinePattern.Solid));
    connect(r1.n, r2.n) annotation(
      Line(points = {{40, 30}, {50, 30}, {50, 10}, {40, 10}}, color = {0, 0, 255}, pattern = LinePattern.Solid));
    connect(r2.n, r3.n) annotation(
      Line(points = {{40, 10}, {50, 10}, {50, -10}, {40, -10}}, color = {0, 0, 255}, pattern = LinePattern.Solid));
    connect(r3.p, line.n[3]) annotation(
      Line(points = {{20, -10}, {0, -10}, {0, 10}}, color = {0, 0, 255}, pattern = LinePattern.Solid));
    connect(v1.n, line.p[1]) annotation(
      Line(points = {{-40, 30}, {-20, 30}, {-20, 10}}, color = {0, 0, 255}, pattern = LinePattern.Solid));
    connect(v2.n, line.p[2]) annotation(
      Line(points = {{-40, 10}, {-20, 10}}, color = {0, 0, 255}, pattern = LinePattern.Solid));
    connect(v3.n, line.p[3]) annotation(
      Line(points = {{-40, -10}, {-20, -10}, {-20, 10}}, color = {0, 0, 255}, pattern = LinePattern.Solid));
    connect(v1.p, rgUL.p) annotation(
      Line(points = {{-60, 30}, {-70, 30}, {-70, -20}}, color = {0, 0, 255}, pattern = LinePattern.Solid));
    connect(v2.p, rgUL.p) annotation(
      Line(points = {{-60, 10}, {-70, 10}, {-70, -20}}, color = {0, 0, 255}, pattern = LinePattern.Solid));
    connect(v3.p, rgUL.p) annotation(
      Line(points = {{-60, -10}, {-70, -10}, {-70, -20}}, color = {0, 0, 255}, pattern = LinePattern.Solid));
    connect(body.p, line.n[4]) annotation(
      Line(points = {{0, -20}, {0, 10}}, color = {0, 0, 255}, pattern = LinePattern.Solid));
    connect(body.n, ground.p) annotation(
      Line(points = {{0, -40}, {0, -50}, {-70, -50}}, color = {0, 0, 255}, pattern = LinePattern.Solid));
    connect(body.n, r3.n) annotation(
      Line(points = {{0, -40}, {0, -50}, {50, -50}, {50, -10}, {40, -10}}, color = {0, 0, 255}, pattern = LinePattern.Solid));
    connect(Conn1.p, line.p[4]) annotation(
      Line(points = {{-20, -20}, {-20, 10}}, color = {0, 0, 255}, pattern = LinePattern.Solid));
    connect(Conn1.n, ground.p) annotation(
      Line(points = {{-20, -40}, {-20, -50}, {-70, -50}}, color = {0, 0, 255}, pattern = LinePattern.Solid));
    annotation(
      Diagram(coordinateSystem(extent = {{-80, -80}, {60, 40}})),
      Documentation(info = "<html>
<p>
This example shows the usage of <a href=\"modelica://Modelica.Electrical.Analog.Lines.M_OLine\">M_OLine</a>
in a overhead PowerLine. The line geometry is the one shown in figure 4.11 of
[<a href=\"modelica://Modelica.Electrical.Analog.UsersGuide.References\">Theory Book</a>].
</p>
<p>
A common rule-of-thumb states that for steady-state analysis at industrial frequency each
segment (PI model) should be not longer than one tenth of the wave length. For 50&nbsp;Hz,
the wavelength is 6000&nbsp;km, therefore each pi should not overcome 60&nbsp;km. This is
why for our 100&nbsp;km line we have chosen just two trunks.
</p>
<p>
This example has some limitations:
</p>
<ul>
  <li>
    the results are valid only for their steady-state trend. If we want to simulate transients
    effectively, the number of trunks must be enlarged consistently. See example
    &quot;CompareLineTrunks&quot;.
  </li>
  <li>
    since it is based on <a href=\"modelica://Modelica.Electrical.Analog.Lines.M_OLine\">M_OLine</a>,
    which requires zero-valued resistances on the off-diagonal elements of Z impedance matrix,
    it does not take into account the influence of the soil resistance into Z matrix; however,
    for simulations like this one in which the current through the soil is negligible and we
    analyse the results at power frequency where soil resistance is very low, this implies
    negligible errors.
  </li>
</ul>
<p>
The considered line is a three-phase line with a fence all along the line (conductor #4).
An interesting result of this example is the fence steady-state voltage due to the capacitive
coupling with the other wires.
</p>

<h5>Test 1:</h5>
<p>
Simulate the model as it is: a person (body) is supposed to be touching the fence, and
absorbs the current <code>body.i</code>, around 1&nbsp;A (RMS), by large lethal.
</p>

<h5>Test 2:</h5>
<p>
Simulate the model changing the body resistance to 1e6&nbsp;&ohm; or removing body: the steady-state
fence voltage <code>line.n[4].v</code> has changed from 1429&nbsp;V to 5594&nbsp;V (peak).
So there exists a marked difference between the open-circuit voltage and actual body
voltage which, although being much lower, is still very dangerous.
</p>

<h5>Test 3:</h5>
<p>
Simulate the original model reducing the line length to 3&nbsp;km. The body current has fallen
to around 30&nbsp;mA (RMS), which is considered the maximum safe current a body can resist over
5&nbsp;s.
</p>
<p>
Note that the fence voltage before the human contact is roughtly the same for all lengths while
the current flowing when contact occurs is very dependent on the line length.
</p>
</html>", revisions = "<html>
<p>June, 2023 Massimo Ceraolo of the University of Pisa </p>
<p>originally created </p>
</html>"),
      experiment(StopTime = 0.04, Interval = 1e-05, Tolerance = 1e-06),
      Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}})));
  end PowerLineWithFenceSC;

  model CompareCmatrix
    extends Modelica.Icons.Example;
    //********************************************

    model MTL_Cmatrix "Y-matrix with capacitive coupling for a multi-conductor line"
      parameter Integer n(final min = 1) = 3 "Number of conductors";
      parameter Types.CapacitancePerUnitLength C[n, n] "Capacitance matrix; off-diagonal elements are negative";
      parameter Modelica.Units.SI.Length len = 100e3 "line length";
      Modelica.Electrical.Analog.Interfaces.PositivePin p[n] "Vector pin" annotation(
        Placement(transformation(extent = {{-110, 0}, {-90, 80}}), iconTransformation(extent = {{-110, 0}, {-90, 80}})));
      Modelica.Electrical.Analog.Interfaces.NegativePin pn "Reference pin" annotation(
        Placement(transformation(extent = {{-112, -84}, {-90, -62}}), iconTransformation(extent = {{-112, -84}, {-90, -62}})));
      Modelica.Units.SI.Voltage v[n](each start = 0, each fixed = true) "Conductor Voltages with respect to reference";
      Modelica.Units.SI.Current i[n] "Current through conductors";
    equation
      for j in 1:n loop
        v[j] = p[j].v - pn.v;
        i[j] = p[j].i;
      end for;
      0 = sum(p[j].i for j in 1:n) + pn.i;
      i = len*C*der(v);
      annotation(
        Icon(graphics = {Rectangle(extent = {{-100, 100}, {100, -100}}, lineColor = {0, 0, 0}, fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid), Text(extent = {{-96, 8}, {98, -14}}, textColor = {0, 0, 0}, textStyle = {TextStyle.Bold}, textString = "i=C*der(v)"), Text(extent = {{-94, 136}, {96, 118}}, textColor = {0, 0, 255}, textString = "%name")}));
    end MTL_Cmatrix;

    //********************************************
    Modelica.Electrical.Analog.Basic.Ground ground annotation(
      Placement(visible = true, transformation(extent = {{90, -80}, {110, -60}}, rotation = 0)));
    parameter Types.ResistancePerUnitLength Rcomp[div(g.n*(g.n + 1), 2)](each fixed = false) "Compact resistance matrix";
    parameter Types.ReactancePerUnitLength Xcomp[div(g.n*(g.n + 1), 2)](each fixed = false) "Compact reactance matrix";
    parameter Types.InductancePerUnitLength Lcomp[div(g.n*(g.n + 1), 2)](each fixed = false) "Compact inductance";
    parameter Types.CapacitancePerUnitLength C[g.n, g.n](each fixed = false);
    parameter Types.CapacitancePerUnitLength Ccomp[div(g.n*(g.n + 1), 2)](each fixed = false);
    parameter Modelica.Units.SI.Length len = 100e3;
    parameter Real CL[div(g.n*(g.n + 1), 2)] = Ccomp*len;
    MTL_Cmatrix MTL_C(n = g.n, C = C, len = len) annotation(
      Placement(transformation(extent = {{100, 0}, {120, 20}})));
    Modelica.Electrical.Analog.Sources.SineVoltage v2(V = v1.V, phase = -2.0943951023932, f = v1.f) annotation(
      Placement(visible = true, transformation(origin = {50, 10}, extent = {{10, 10}, {-10, -10}}, rotation = 180)));
    Modelica.Electrical.Analog.Sources.SineVoltage v3(V = v1.V, phase = 2.0943951023932, f = v1.f) annotation(
      Placement(visible = true, transformation(origin = {50, -20}, extent = {{10, 10}, {-10, -10}}, rotation = 180)));
    Modelica.Electrical.Analog.Sources.SineVoltage v1(V = 345e3*sqrt(2/3), f = 60) annotation(
      Placement(visible = true, transformation(origin = {50, 40}, extent = {{10, 10}, {-10, -10}}, rotation = 180)));
    Modelica.Electrical.Analog.Basic.Resistor r1(R = r1n.R) annotation(
      Placement(visible = true, transformation(origin = {80, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Electrical.Analog.Basic.Resistor r2(R = r1.R) annotation(
      Placement(visible = true, transformation(origin = {80, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Electrical.Analog.Basic.Resistor r3(R = r1.R) annotation(
      Placement(visible = true, transformation(origin = {80, -20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Electrical.Analog.Basic.Capacitor c12(C = CL[2]) annotation(
      Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = -90, origin = {-10, 30})));
    Modelica.Electrical.Analog.Basic.Capacitor c23(v(fixed = false), C = CL[5]) annotation(
      Placement(transformation(extent = {{10, 10}, {-10, -10}}, rotation = 90, origin = {-10, -10})));
    Modelica.Electrical.Analog.Basic.Capacitor c10(v(fixed = true), C = CL[1]) annotation(
      Placement(transformation(extent = {{10, 10}, {-10, -10}}, rotation = 90, origin = {-60, -40})));
    Modelica.Electrical.Analog.Basic.Capacitor c13(v(fixed = false), C = CL[3]) annotation(
      Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = -90, origin = {10, 10})));
    Modelica.Electrical.Analog.Basic.Capacitor c30(v(fixed = true), C = CL[6]) annotation(
      Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = -90, origin = {-20, -40})));
    Modelica.Electrical.Analog.Basic.Capacitor c20(v(fixed = true), C = CL[4]) annotation(
      Placement(transformation(extent = {{-10, -10}, {10, 10}}, rotation = -90, origin = {-40, -40})));
    Modelica.Electrical.Analog.Basic.Ground ground1 annotation(
      Placement(visible = true, transformation(extent = {{-30, -90}, {-10, -70}}, rotation = 0)));
    Modelica.Electrical.Analog.Sources.SineVoltage v2a(V = v1.V, phase = -2.0943951023932, f = v1.f) annotation(
      Placement(visible = true, transformation(origin = {-110, 10}, extent = {{10, 10}, {-10, -10}}, rotation = 180)));
    Modelica.Electrical.Analog.Sources.SineVoltage v3a(V = v1.V, phase = 2.0943951023932, f = v1.f) annotation(
      Placement(visible = true, transformation(origin = {-110, -20}, extent = {{10, 10}, {-10, -10}}, rotation = 180)));
    Modelica.Electrical.Analog.Sources.SineVoltage v1a(V = 345e3*sqrt(2/3), f = 60) annotation(
      Placement(visible = true, transformation(origin = {-110, 40}, extent = {{10, 10}, {-10, -10}}, rotation = 180)));
    parameter PowerLineMC.Functions.LineGeometry g(n = 3, x = {0, -3.048, 3.048}, y = {12.192, 12.192, 12.192}, r = 1e-3/2*{12.7, 12.7, 12.7}, R1 = 1e-3*{0.348, 0.348, 0.348}, k_s = {0.287, 0.287, 0.287}, f = 60) annotation(
      Placement(visible = true, transformation(origin = {-60, -128}, extent = {{80, 60}, {100, 80}}, rotation = 0)));
    Modelica.Electrical.Analog.Basic.Resistor r1n(R = 1000.0) annotation(
      Placement(visible = true, transformation(origin = {-80, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Electrical.Analog.Basic.Resistor r2n(R = r1n.R) annotation(
      Placement(visible = true, transformation(origin = {-80, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Electrical.Analog.Basic.Resistor r3n(R = r1n.R) annotation(
      Placement(visible = true, transformation(origin = {-80, -20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  initial algorithm
    (Ccomp, C) := Functions.lineCmatrix(n = g.n, x = g.x, y = g.y, r = g.r);
    (Rcomp, Xcomp, Lcomp) := Functions.lineZmatrix(n = g.n, x = g.x, y = g.y, r = g.r, R1 = g.R1, k_s = g.k_s, rho = g.rho, f = g.f);
  equation
    connect(v1.n, r1.p) annotation(
      Line(points = {{60, 40}, {70, 40}}, color = {0, 0, 255}));
    connect(v2.n, r2.p) annotation(
      Line(points = {{60, 10}, {70, 10}}, color = {0, 0, 255}));
    connect(v3.n, r3.p) annotation(
      Line(points = {{60, -20}, {70, -20}}, color = {0, 0, 255}));
    connect(r1.n, MTL_C.p[1]) annotation(
      Line(points = {{90, 40}, {96, 40}, {96, 18}, {100, 18}, {100, 14}}, color = {0, 0, 255}));
    connect(r2.n, MTL_C.p[2]) annotation(
      Line(points = {{90, 10}, {90, 14}, {100, 14}}, color = {0, 0, 255}));
    connect(r3.n, MTL_C.p[3]) annotation(
      Line(points = {{90, -20}, {96, -20}, {96, 10}, {100, 10}, {100, 14}}, color = {0, 0, 255}));
    connect(v1.p, v2.p) annotation(
      Line(points = {{40, 40}, {30, 40}, {30, 10}, {40, 10}}, color = {0, 0, 255}));
    connect(v2.p, v3.p) annotation(
      Line(points = {{40, 10}, {30, 10}, {30, -20}, {40, -20}}, color = {0, 0, 255}));
    connect(MTL_C.pn, ground.p) annotation(
      Line(points = {{99.9, 2.7}, {100, 2.7}, {100, -60}}, color = {0, 0, 255}));
    connect(c30.n, ground1.p) annotation(
      Line(points = {{-20, -50}, {-20, -70}}, color = {0, 0, 255}));
    connect(c20.n, ground1.p) annotation(
      Line(points = {{-40, -50}, {-40, -60}, {-20, -60}, {-20, -70}}, color = {0, 0, 255}));
    connect(c10.n, ground1.p) annotation(
      Line(points = {{-60, -50}, {-60, -60}, {-20, -60}, {-20, -70}}, color = {0, 0, 255}));
    connect(c10.p, c12.p) annotation(
      Line(points = {{-60, -30}, {-60, 40}, {-10, 40}}, color = {0, 0, 255}));
    connect(c12.n, c23.p) annotation(
      Line(points = {{-10, 20}, {-10, 0}}, color = {0, 0, 255}));
    connect(v1a.p, v2a.p) annotation(
      Line(points = {{-120, 40}, {-120, 10}}, color = {0, 0, 255}));
    connect(v2a.p, v3a.p) annotation(
      Line(points = {{-120, 10}, {-120, -20}}, color = {0, 0, 255}));
    connect(c12.p, c13.p) annotation(
      Line(points = {{-10, 40}, {10, 40}, {10, 20}}, color = {0, 0, 255}));
    connect(c13.n, c23.n) annotation(
      Line(points = {{10, 0}, {10, -20}, {-10, -20}}, color = {0, 0, 255}));
    connect(v1a.n, r1n.p) annotation(
      Line(points = {{-100, 40}, {-90, 40}}, color = {0, 0, 255}));
    connect(r1n.n, c12.p) annotation(
      Line(points = {{-70, 40}, {-10, 40}}, color = {0, 0, 255}));
    connect(v3a.n, r3n.p) annotation(
      Line(points = {{-100, -20}, {-90, -20}}, color = {0, 0, 255}));
    connect(v2a.n, r2n.p) annotation(
      Line(points = {{-100, 10}, {-90, 10}}, color = {0, 0, 255}));
    connect(r2n.n, c23.p) annotation(
      Line(points = {{-70, 10}, {-10, 10}, {-10, 0}}, color = {0, 0, 255}));
    connect(r3n.n, c23.n) annotation(
      Line(points = {{-70, -20}, {-10, -20}}, color = {0, 0, 255}));
    connect(ground1.p, v3a.p) annotation(
      Line(points = {{-20, -70}, {-20, -60}, {-120, -60}, {-120, -20}}, color = {0, 0, 255}));
    connect(ground.p, v3.p) annotation(
      Line(points = {{100, -60}, {100, -40}, {30, -40}, {30, -20}, {40, -20}}, color = {0, 0, 255}));
    connect(c20.p, c23.p) annotation(
      Line(points = {{-40, -30}, {-40, 10}, {-10, 10}, {-10, 0}}, color = {0, 0, 255}));
    connect(c30.p, c23.n) annotation(
      Line(points = {{-20, -30}, {-20, -20}, {-10, -20}}, color = {0, 0, 255}));
    annotation(
      Diagram(coordinateSystem(extent = {{-120, 60}, {120, -100}})),
      Documentation(info = "<html>
<p>
MSL uses, inside <a href=\"modelica://Modelica.Electrical.Analog.Lines.M_OLine\">M_OLine</a>,
a network of physical capacitors to model capacitive effects across conductors.
</p>
<p>
Instead, in typical PowerSystems approach, a capacitance matrix is used, corresponding to
the vector equation:
</p>
<p>
<strong>i</strong>=<em><strong>C&nbsp;</strong></em>d<strong>v</strong>/dt&nbsp;
</p>
<p>
where
</p>
<ul>
  <li>
    <em><strong>i</strong></em> is the vector of currents entering a multiple-conductor lines
    (n conductors with a return wire gives rise to n currents),
  </li>
  <li>
    <em><strong>v</strong></em> is the vector of voltages at a port of a multiple-conductor
    lines (n conductors with a return wire gives rise to n voltages).
  </li>
</ul>
<p>
The <em><strong>C</strong></em> matrix is symmetrical, and therefore it has n(n+1) independent
terms. the diagonal terms of this matrix are negative.
</p>
<p>
The capacitive coupling across conductors can be modelled equivalently through a network of
n(n+1) capacitors, and simple conversion formulas exist, which can for instance be found in
[<a href=\"modelica://Modelica.Electrical.Analog.UsersGuide.References\">Theory Book</a>],
page 3-12.
</p>
<p>
The support function
<a href=\"modelica://Modelica.Electrical.Analog.Lines.Functions.lineCmatrix\">LineCmatrix</a>
already contains these formulas and puts these capacitance values in an array called
<strong>Ccompact</strong>.
</p>
<p>
This example is to show the perfect equivalence of capacitor matrix and Ccompact values
(here stored in an array named <strong>Ccomp</strong>). Just check for equality
<code>r1n.i</code> and <code>r1.i</code>.
</p>
</html>", revisions = "<html>
<p><em>June, 2023 Massimo Ceraolo of the University of Pisa</em></p>
<p><em>originally created</em></p>
</html>"),
      experiment(StopTime = 0.04, Interval = 1e-05, Tolerance = 1e-06),
      Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}})));
  end CompareCmatrix;

  package Functions "Functions and records"
    extends Modelica.Icons.FunctionsPackage;

    record LineGeometry
      extends Modelica.Icons.Record;
      parameter Integer n = 3 "Number of wires" annotation(
        Evaluate = true);
      parameter Modelica.Units.SI.Length x[n] "Horizontal abscissas of conductors";
      parameter Modelica.Units.SI.Length y[n] "Vertical abscissas of conductors";
      parameter Modelica.Units.SI.Length r[n] "Conductor radii";
      parameter Types.ResistancePerUnitLength R1[n] "Resistance per length of conductors";
      parameter Real k_s[n] = fill(0.7, n) "Ratio of GMR to actual conductor radius";
      parameter Modelica.Units.SI.Resistivity rho = 100 "Earth resistivity";
      parameter Modelica.Units.SI.Frequency f = 50 "Line frequency";
      annotation(
        Documentation(info = "<html>
<p>This record contains the line geometry and resistances of a multi-conductor line.</p>
<p>It will be used by functions lineZmatrix and lineCmatrix that compute longitudinal (flow) and transverse (cross) line matrices.</p>
<p>Even though frequency is not part of the line geometry, it is a parameter needed by lineZmatrix and therefore it is included here.</p>
</html>", revisions = "<html>
<p><em>July, 2023</em></p>
<p>Original implementation by Massimo Ceraolo of the University of Pisa </p>
</html>"));
    end LineGeometry;

    function lineCmatrix "Compute matrix of transverse capacitances per metre of a multi-conductor line"
      extends Modelica.Icons.Function;
      import Modelica.Constants.epsilon_0;
      import Modelica.Constants.pi;
      input Integer n "Number of conductors";
      input Modelica.Units.SI.Length x[n] "Horizontal abscissas of conductors";
      input Modelica.Units.SI.Length y[n] "Vertical abscissas of conductors";
      input Modelica.Units.SI.Radius r[n] "Conductors radii (m)";
      //  output Real Cflat[n, n](each final unit="F/m") "Matrix with capacitances of a network reproducing the behaviour of C";
      output Types.CapacitancePerUnitLength Ccompact[div(n*(n + 1), 2)] "Vector of capacitances of network of capacitors equivalent to C";
      output Types.CapacitancePerUnitLength C[n, n] "Capacitance matrix with negative off-diagonal conductances";
    protected
      constant Real K(final unit = "m/F") = 1/(2*pi*epsilon_0);
      Real p[n, n](each final unit = "m/F") "Maxwell's potential matrix";
      Types.CapacitancePerUnitLength Cflat[n, n] "Matrix with capacitances of a network reproducing the behaviour of C";
      Modelica.Units.SI.Distance D "Generic larger distance";
      Modelica.Units.SI.Distance d "Generic smaller distance";
      Integer k;
      String sC;
    algorithm
// Diagonal elements of the potential matrix:
      for i in 1:n loop
        p[i, i] := K*log(2*y[i]/r[i]);
      end for;
// Out-of diagonal elements of the potential matrix:
      for i in 1:n loop
        for jj in 1:i - 1 loop
          d := sqrt((x[i] - x[jj])^2 + (y[i] - y[jj])^2);
          D := sqrt((x[i] - x[jj])^2 + (y[i] + y[jj])^2);
          p[i, jj] := K*log(D/d);
        end for;
      end for;
      for i in 1:n loop
        for jj in i + 1:n loop
          p[i, jj] := p[jj, i];
        end for;
      end for;
// The capacitance function is the inverse of the matrix of potentials
      C := Modelica.Math.Matrices.inv(p);
// The off diagonal elements of Cflat are the opposite of corresponding elements of C.
// We first generate the opposite of C, then will redefine its diagonal elements:
      Cflat := -C;
// The diagonal elements of Cflat contain the sum of all the values the corresponding row in C:
      for i in 1:n loop
        Cflat[i, i] := sum(C[i, j] for j in 1:n);
      end for;
// Select the elements needed by M_Oline in a vector which it can use directly
      k := 0;
      for i in 1:n loop
        for j in i:n loop
          k := k + 1;
          Ccompact[k] := Cflat[i, j];
        end for;
      end for;
      annotation(
        Documentation(info = "<html>
<p>This function computes Capacitances of multi-conductor transmission lines, according to the formulas as reported in [<a href=\"modelica://Modelica.Electrical.Analog.UsersGuide.References\">Cerolo2018</a>, Appendix]. The results obtained with this function have been checked with Fig. 4.1 of [<a href=\"modelica://Modelica.Electrical.Analog.UsersGuide.References\">Theory Book</a>], with good compliance. </p>
<p>Internally, it computes the <strong>C</strong> matrix, which corresponds to the formulas <strong>V</strong>&nbsp;=&nbsp;<strong>YI</strong>, <strong>Y</strong>&nbsp;=&nbsp;&omega;<strong>C</strong> where </p>
<ul>
<li><strong>V</strong> is the vector of voltages between conductors and the reference (the return conductor, usually ground), </li>
<li><strong>I</strong>&nbsp;is the vector of transverse currents (between conductors and the return conductor, usually ground) due to the capacitive coupling between the conductors, </li>
<li><strong>Y</strong> is the matrix of transverse admittances of the multiconductor line (S/m), </li>
<li>&omega; is the angular frequency when constant-frequency steady-state operation of the line is considered </li>
</ul>
<p>This matrix&nbsp;<strong>C</strong>, has always negative off-diagonal values, and positive diagonal values. </p>
<p>From <strong>C</strong> matrix, the internal <strong>Cflat</strong> matrix is computed, containing physical capacitors that can be imagined between conductors to model capacitive effects. For instance C12 is the capacitance (per unit length) to be put between conducturors 1 and 2. The output array <strong>Ccompact</strong> contains the elements of the <strong>Cflat</strong> matrix ordered as described in the <a href=\"modelica://Modelica.Electrical.Analog.Lines.M_OLine\">M_OLine</a> model, and is used in example <a href=\"modelica://Modelica.Electrical.Analog.Examples.Lines.PowerLineWithFence\">Examples.Lines.PowerLineWithFence</a> in conjunction with M_OLine. </p>
<p>For an example on how to use this function, consider model Electrical.Analog.Examples.Lines.CompareCmatrix. </p>
</html>", revisions = "<html>
<p><em>July, 2023</em> </p>
<p>Original implementation by Massimo Ceraolo of the University of Pisa </p>
</html>"));
    end lineCmatrix;

    function lineZmatrix "Compute matrix of longitudinal impedances per metre of a multi-conductor line"
      extends Modelica.Icons.Function;
      import Modelica.ComplexMath;
      import Modelica.Utilities.Streams;
      import Modelica.Constants.pi;
      import Modelica.Constants.mu_0;
      import Modelica.ComplexMath.j;
      input Integer n "Number of conductors in line";
      input Modelica.Units.SI.Length x[n] "Horizontal abscissas of conductors";
      input Modelica.Units.SI.Length y[n] "Vertical abscissas of conductors";
      input Modelica.Units.SI.Length r[n] "Conductors radii";
      input Types.ResistancePerUnitLength R1[n] "Conductors lineic resistance";
      input Real k_s[n] = fill(0.7, n) "Ratio of equivalent shell radius to actual radius";
      //  in case of cylindric conductor this is equal to exp(-mu_r/4)=exp(-0.25)
      input Modelica.Units.SI.Resistivity rho = 100 "Ground resistivity";
      input Modelica.Units.SI.Frequency f = 50 "Frequency";
      output Types.ResistancePerUnitLength Rcomp[div(n*(n + 1), 2)] "Compact resistance matrix";
      output Types.ReactancePerUnitLength Xcomp[div(n*(n + 1), 2)] "Compact reactance matrix";
      output Types.InductancePerUnitLength Lcomp[div(n*(n + 1), 2)] "Compact inductance)";
    protected
      Modelica.Units.SI.Distance D "Generic larger distance";
      Modelica.Units.SI.Distance d "Generic smaller distance";
      Real Theta "Theta angle for Carson's formulas";
      Real L "generic inductance";
      Real P "P coefficient from Carson's original paper";
      Real Q "Q coefficient from Carson's original paper";
      Complex Z[n, n] "Computed matrix (ohm/m)";
      final parameter Real w = 2*pi*f;
      parameter Real a0 = sqrt(w*mu_0/rho) "multiplies D in a-parameter from EMTPs Theory Book";
      // a0 is written, equivalently, as 4*pi*sqrt(5)*1e-4*sqrt(f/rho) in EMTP Theory book
      //(reference available in Electrical.Analog.References)
      //
      //  a0*D=a; a is the same as "r" in Carson's paper
      // where D=2h for self impedance, D=Dik for mutual impedance
      Real a;
      Integer k;
    algorithm
//self impedances with ground return:
      for i in 1:n loop
        a := a0*2*y[i];
        assert(a < 0.25, "parameter a=" + String(a) + " >0.25 will result in inadequate precision of Carson's formulas", AssertionLevel.warning);
        P := pi/8 - sqrt(2)/6*a + a^2/16*(0.672784 + log(2/a));
        Q := (-0.03860784) + 0.5*log(2/a) + sqrt(2)/6*a - pi/64*a*a;
        L := mu_0/(2*pi)*log(2*y[i]/(k_s[i]*r[i])) "L when soil is pure conductive";
        Z[i, i] := R1[i] + 4*w*P/1.e7 + j*w*(L + 4*Q/1.e7);
//Carson (25) and (30)
      end for;
//mutual impedances with ground return:
      for i in 1:n loop
        for jj in 1:i - 1 loop
          Theta := atan(abs((x[i] - x[jj])/(y[i] + y[jj])));
          d := sqrt((x[i] - x[jj])^2 + (y[i] - y[jj])^2);
//distance between conductors i and jj
          D := sqrt((x[i] - x[jj])^2 + (y[i] + y[jj])^2);
//distance between conductor i and the image of conductor jj
          a := a0*D;
          P := pi/8 - sqrt(2)/6*a*cos(Theta) + a^2/16*(0.672784 + log(2/a))*cos(2*Theta) + a^2/16*Theta*sin(2*Theta);
          Q := (-0.03860784) + 0.5*log(2/a) + sqrt(2)/6*a*cos(Theta) - pi/64*a*a*cos(2*Theta);
          L := mu_0/(2*pi)*log(D/d) "L when soil is pure conductive";
          Z[i, jj] := 4*w*P/1.e7 + j*w*(L + 4*Q/1.e7);
//Carson (26) and (31)
        end for;
      end for;
      k := 0;
      for i in 1:n loop
        for jj in i:n loop
          k := k + 1;
          Rcomp[k] := Z[jj, i].re;
          Xcomp[k] := Z[jj, i].im;
          Lcomp[k] := Xcomp[k]/w;
        end for;
      end for;
      annotation(
        Documentation(info = "<html>
<p>This function computes resistances, reactances, inductances of multi-conductor transmission lines, taking into account ground characteristics according to John Carson&apos;s formulas as reported in [<a href=\"modelica://Modelica.Electrical.Analog.UsersGuide.References\">Cerolo2018</a>, Appendix], but upgraded in accuracy (more digits in constants and more added terms in the evaluation of Q). </p>
<p>When used in conjunction with <a href=\"modelica://Modelica.Electrical.Analog.Lines.M_OLine\">M_OLine</a> from its output just the inductances are used; the resistances could instead be used when ground resistance is to be taken into account. </p>
<p>The results obtained with this function have been checked with Fig. 4.11 of [<a href=\"modelica://Modelica.Electrical.Analog.UsersGuide.References\">Theory Book</a>], in model <a href=\"modelica://Modelica.Electrical.Analog.Examples.Lines.PowerLineWithFence\">PowerLineWithFence</a>, with good agreement. </p>
<p>The output arrays contain the elements of the Z matrix ordered as described in the <a href=\"modelica://Modelica.Electrical.Analog.Lines.M_OLine\">M_OLine</a> model and are used in example Examples.Lines.PowerLineWithFence in conjunction with M_OLine. </p>
<p>Parameter <span style=\"font-family: Courier New;\">k_s</span> takes into account the conductor&apos;s skin effect; it is 0.778 for a solid non-magnetic conductor, can be between 0.35 and 0.8 for real-life power line conductors (they are stranded and often have an iron core). </p>
<p>According to Carson&apos;s theory the line impedances depend on the frequency of signal, and therefore this function is run for the input given frequency. This is mainly because <span style=\"font-family: MS Shell Dlg 2;\">there is an important skin effect of the ground that causes ground return resistances and inductances to respectively increase and decrease as frequency increases. Nevertheless, values such as those obtained from these formulas are commonly used in power system transients, the steady-state frequency (e.g. 50 or 60 Hz) being used to compute Z matrix from this function.</span></p>
<p>The formulas inside have adequate precision only if the requested computation frequency is within given limits (up to several hundred Hertz); if a larger than acceptable frequency is requested, given the line geometry, a warning is issued. </p>
</html>", revisions = "<html>
<p><em>July, 2023</em> </p>
<p>Original implementation by Massimo Ceraolo of the University of Pisa </p>
</html>"));
    end lineZmatrix;

    model TestLineZmatrix
      extends Modelica.Icons.Example;
      import Modelica.Utilities.Streams.print;
      parameter LineGeometry g(n = 4, x = {0, -3.048, 3.048, -9.144}, y = {12.192, 12.192, 12.192, 3.048}, r = 1e-3/2*{12.7, 12.7, 12.7, 4.064}, R1 = 1e-3*{0.348, 0.348, 0.348, 1.802}, k_s = {0.287, 0.287, 0.287, 0.779}, f = 60) annotation(
        Placement(transformation(extent = {{-10, -6}, {10, 14}})));
      String sC;
      Types.ResistancePerUnitLength Rcomp[div(g.n*(g.n + 1), 2)] "Compacted matrix";
      Types.ReactancePerUnitLength Xcomp[div(g.n*(g.n + 1), 2)] "Compacted matrix";
    protected
      Integer k;
    algorithm
      when initial() then
        (Rcomp, Xcomp) := lineZmatrix(n = g.n, x = g.x, y = g.y, r = g.r, R1 = g.R1, k_s = g.k_s, rho = g.rho, f = g.f);
        print("\n *****              Using lineZmatrix, RESULTS in ohm/km              *****");
        print(" *** (one row per matrix row; numbers should be intended right-aligned) ***");
        k := 0;
        for i in 1:g.n loop
//matrix row
          sC := "";
          for j in 1:g.n - i + 1 loop
// matrix column
            k := k + 1;
            sC := sC + String(1000*Rcomp[k]) + "+j" + String(1000*Xcomp[k]) + "   ";
          end for;
          print(sC);
        end for;
      end when;
      annotation(
        Documentation(info = "<html>
<p>This model tests the Z matrix as computed with function lineZmatrix, with the geometry of fig. 4.11 of [<a href=\"modelica://Modelica.Electrical.Analog.UsersGuide.References\">Theory Book</a>]. </p>
<p>The results are given textually in the log and show a good agreement with the results shown in [<a href=\"modelica://Modelica.Electrical.Analog.UsersGuide.References\">Theory Book</a>].</p>
<p>Expected output:</p>
<pre>
*****                Using lineZmatrix, RESULTS in ohm/km                *****
***   (one row per matrix row; numbers should be intended right-aligned)   ***
0.405445+j0.986077 0.0574443+j0.426468 0.0574443+j0.426468 0.058076+j0.316811 
0.405445+j0.986077 0.0574408+j0.374207 0.0580827+j0.329078 
0.405445+j0.986077 0.0580667+j0.304429 
1.86076+j0.995306 
... &quot;TestLineZmatrix.mat&quot; creating (simulation result file)
</pre>
</html>", revisions = "<html>
<p><em>July, 2023</em> </p>
<p>Original implementation by Massimo Ceraolo of the University of Pisa </p>
</html>"),
        experiment(StopTime = 0, __Dymola_Algorithm = "Dassl"));
    end TestLineZmatrix;

    model TestLineCmatrix
      extends Modelica.Icons.Example;
      import Modelica.Utilities.Streams.print;
      parameter LineGeometry g(n = 4, x = {0, -3.048, 3.048, -9.144}, y = {12.192, 12.192, 12.192, 3.048}, r = 1e-3/2*{12.7, 12.7, 12.7, 4.064}, R1 = 1e-3*{0.348, 0.348, 0.348, 1.802}) annotation(
        Placement(transformation(extent = {{-10, -6}, {10, 14}})));
      String sC;
      Types.CapacitancePerUnitLength Ccomp[div(g.n*(g.n + 1), 2)] "Compacted matrix";
      Types.CapacitancePerUnitLength C[g.n, g.n] "Full C matrix";
    protected
      Integer k;
    algorithm
      when initial() then
        (Ccomp, C) := lineCmatrix(n = g.n, x = g.x, y = g.y, r = g.r);
        print("\n ***** Full C matrix from  lineCmatrix in nF/km *****");
        print(" ***       (only half: matrix is symmetric)       ***");
        for i in 1:g.n loop
//matrix row
          sC := "";
          for j in 1:i loop
// matrix column
            sC := sC + String(1e12*C[i, j]) + "  ";
          end for;
          print(sC);
        end for;
        print("\n *****        matrix of capacitor capacitances from lineCmatrix in nF/km         *****");
        print(" ***  (Cij is capacitance of capacitor between i and j to mime C-matrix behaviour  ***");
        k := 0;
        for i in 1:g.n loop
//matrix row
          sC := "";
          for j in 1:g.n - i + 1 loop
// matrix column
            k := k + 1;
            sC := sC + String(1e12*Ccomp[k]) + "  ";
          end for;
          print(sC);
        end for;
      end when;
      annotation(
        Documentation(info = "<html>
<p>This model tests the C matrix as computed with function lineCmatrix, with the geometry of fig. 4.11 of [<a href=\"modelica://Modelica.Electrical.Analog.UsersGuide.References\">Theory Book</a>]. </p>
<p>The results are given textually in the log and show a good agreement with the results shown in [<a href=\"modelica://Modelica.Electrical.Analog.UsersGuide.References\">Theory Book</a>].</p>
<p>The function computes both the C-matrix and the array of physical capacitances miming C-matrix behaviour (see documentation of lineCmatrix for details)</p>
<p>Expected output:</p>
<pre>
***** Full C matrix from lineCmatrix in nF/km *****
***     (only half: matrix is symmetric)        ***
7.57087
-1.62658 7.30876
-1.63038 -0.834876 7.29987
-0.168826 -0.275822 -0.118934 6.97273

*****       matrix of capacitor capacitances from lineCmatrix in nF/km        *****
*** (Cij is capacitance of capacitor between i and j to mime C-matrix behaviour ***
4.14508 1.62658 1.63038 0.168826
4.57148 0.834876 0.275822
4.71569 0.118934
6.40915
</pre>
</html>", revisions = "<html>
<p><em>July, 2023</em> </p>
<p>Original implementation by Massimo Ceraolo of the University of Pisa </p>
</html>"),
        experiment(StopTime = 0, __Dymola_Algorithm = "Dassl"));
    end TestLineCmatrix;
    annotation(
      Documentation(info = "<html>
<p>This package class contains functions lineCmatrix and lineZmatrix that are able to compute C and Z matrices for a multple conductor transmission line (typically a power line), using one of the possible implementations of Carson&apos;s theory; see their documentation for details. </p>
<p>They are used in Modelica.Electrical.Analog.Examples.Lines models CompareCmatrix and PowerLineWIthFence (see their documentation for details). The latter model is an example of usage of Modelica.Electrical.Analog.Lines.M_Oline model.</p>
<p>These functions make usage of the lineGeometry record, which defines the geometry of a multi-conductor transmission line, with possibile ground return.</p>
</html>"));
  end Functions;

  package Types
    type CapacitancePerUnitLength = Real(final quantity = "CapacitancePerUnitLength", final unit = "F/m") "Capacitance per unit length of wire/cable/line";
    type ConductancePerUnitLength = Real(final quantity = "ConductancePerUnitLength", final unit = "S/m") "Conductance per unit length of wire/cable/line";
    type ImpedancePerUnitLength = ResistancePerUnitLength "Impedance per unit length of wire/cable/line";
    type InductancePerUnitLength = Real(final quantity = "InductancePerUnitLength", final unit = "H/m") "Inductance per unit length of wire/cable/line";
    type ReactancePerUnitLength = ResistancePerUnitLength "Reactance per unit length of wire/cable/line";
    type ResistancePerUnitLength = Real(final quantity = "ResistancePerUnitLength", final unit = "Ohm/m") "Resistance per unit length of wire/cable/line";
  end Types;
  annotation(
    uses(Modelica(version = "4.0.0"), Complex(version = "4.0.0")));
end PowerLineMC;
