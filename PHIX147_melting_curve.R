library(ggplot2)
df <- read.csv("/Users/jaren/Desktop/BTEC 640/Bob Horton Project/Experimental_Melting_Curve_PHIX174_normalized_from_paper_and_digitalized.csv")


ggplot(df, aes(x = Temperature, y = RawValue)) +
  geom_line(color = "black", linewidth = 1.2) +
  geom_point(size = 0.7, color = "black") +
  labs(
    title = " \u03C6X174 Experimental Melting Curve",
    x = "T",
    y = expression(dN[c]/dT ~ 'x 10'^-3)
  ) +
  theme_minimal(base_size = 14)



#for decipher plot 
df <- read.csv("/Users/jaren/Desktop/BTEC 640/Bob Horton Project/DECIPHER_phiX174_melt_curve.csv")
ggplot(df, aes(x = Temperature, y = Derivative)) +
  geom_line(color = "black", linewidth = 1.2) +
  geom_point(size = 0.7, color = "black") +
  labs(
    title = " \u03C6X174 Experimental Melting Curve",
    x = "T",
    y = expression(dN[c]/dT ~ 'x 10'^-3)
  ) +
  theme_minimal(base_size = 14)

