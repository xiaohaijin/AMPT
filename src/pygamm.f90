Function pygamm(x)
  Dimension b(8)
  Data b/ -0.57719165, 0.98820589, -0.89705694, 0.91820686, -0.75670408, 0.48219939, -0.19352782, 0.03586834/
  nx = int(x)
  dx = x - nx
  pygamm = 1.
  Do i = 1, 8
    pygamm = pygamm + b(i)*dx**i
  End Do
  If (x<1.) Then
    pygamm = pygamm/x
  Else
    Do ix = 1, nx - 1
      pygamm = (x-ix)*pygamm
    End Do
  End If
  Return
End Function pygamm
