within ;
package BodylightExtension
  package Types
    type MolarReactionRate = Real(
      quantity="MolarReactionRate",
      unit = "mol/(m3.s)",
      displayUnit = "mmol/(l.s)"
      );
    type MolarFluxPerArea = Real(
      quantity    = "MolarFluxPerArea",
      unit        = "mol/(m2.s)",
      displayUnit = "mol/(cm2.s)"
      );
    type SurfaceConcentration = Real(
      quantity="SurfaceConcentration",
      unit="mol/m2",
      displayUnit="mol/cm2"
      );
    type DiffusionCoefficient= Real(
      quantity    = "DiffusionCoefficient",
      unit        = "m2/s",
      displayUnit = "cm2/s"
      );
    type ReactionRateFirstOrder = Real(
      quantity    = "ReactionRateFirstOrder",
      unit        = "1/s"
      );
    type ReactionRateSecondOrder = Real(
      quantity    = "ReactionRateSecondOrder",
      unit        = "m3/(mol.s)",
      displayUnit = "l/(mol.s)"
      );
    type ReactionRateThirdOrder = Real(
      quantity    = "ReactionRateThirdOrder",
      unit        = "m6/(mol2.s)",
      displayUnit = "l2/(mol2.s)"
      );
    package RealIO
      connector MolarReactionRateInput = input MolarReactionRate
      "input MolarReactionRate as connector"
      annotation (defaultComponentName="molarReactionRate",
        Icon(graphics={Polygon(
                points={{-100,100},{100,0},{-100,-100},{-100,100}},
                lineColor={0,0,127},
                fillColor={0,0,127},
                fillPattern=FillPattern.Solid)},
             coordinateSystem(extent={{-100,-100},{100,100}}, preserveAspectRatio=true, initialScale=0.2)),
        Diagram(coordinateSystem(
              preserveAspectRatio=true, initialScale=0.2,
              extent={{-100,-100},{100,100}},
              grid={1,1}), graphics={Polygon(
                points={{0,50},{100,0},{0,-50},{0,50}},
                lineColor={0,0,127},
                fillColor={0,0,127},
                fillPattern=FillPattern.Solid), Text(
                extent={{-10,85},{-10,60}},
                lineColor={0,0,127},
                textString="%name")}),
          Documentation(info="<html>
    <p>
    Connector with one input signal of type MolarFlowRate.
    </p>
    </html>"));
      connector MolarReactionRateOutput= output MolarReactionRate
        "output MolarReactionRate as connector"
      annotation (defaultComponentName="molarReactionRate",
        Icon(coordinateSystem(
              preserveAspectRatio=true,
              extent={{-100,-100},{100,100}},
              grid={1,1}), graphics={Polygon(
                points={{-100,100},{100,0},{-100,-100},{-100,100}},
                lineColor={0,0,127},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid)}),
        Diagram(coordinateSystem(
              preserveAspectRatio=true,
              extent={{-100,-100},{100,100}},
              grid={1,1}), graphics={Polygon(
                points={{-100,50},{0,0},{-100,-50},{-100,50}},
                lineColor={0,0,127},
                fillColor={255,255,255},
                fillPattern=FillPattern.Solid), Text(
                extent={{30,110},{30,60}},
                lineColor={0,0,127},
                textString="%name")}),
          Documentation(info="<html>
  <p>
  Connector with one output signal of type Real.
  </p>
  </html>"));
    end RealIO;

    package Constants
      extends Modelica.Icons.SourcesPackage;


    block MolarFlowRateConst "Constant signal of type MolarFlowRate"
     parameter BodylightExtension.Types.MolarReactionRate k
          "Constant MolarFlowRate output value";
          BodylightExtension.Types.RealIO.MolarReactionRateOutput y "MolarFlowRate constant"
          annotation (Placement(transformation(extent={{40,-10},{60,10}}),
              iconTransformation(extent={{40,-10},{60,10}})));
    equation
          y=k;
      annotation (defaultComponentName="molarFlowRate",
                 Diagram(coordinateSystem(extent={{-40,-40},{40,40}})), Icon(
            coordinateSystem(extent={{-40,-40},{40,40}}, preserveAspectRatio=false),
                graphics={
            Rectangle(extent={{-40,40},{40,-40}},
              lineColor={0,0,0},
                  radius=10,
              fillColor={236,236,236},
                              fillPattern=FillPattern.Solid),
            Text( extent={{-100,-44},{100,-64}},
              lineColor={0,0,0},
                      fillColor={236,236,236},
              fillPattern=FillPattern.Solid,
                  textString="%name"),
            Text(         extent={{-40,10},{40,-10}},
              lineColor={0,0,0},
                  fillColor={236,236,236},
              fillPattern=FillPattern.Solid,
                      textString="Const")}));
    end MolarFlowRateConst;
    end Constants;

  end Types;
  annotation (uses(Bodylight(version="1.0")));
end BodylightExtension;
