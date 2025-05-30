#include <TColor.h>

/// Color schemes, trying to be as color-vision-deficiency friendly as possible (hexadecimal codes)

// new official CMS recommendations, see CAT group's doc https://cms-analysis.docs.cern.ch/guidelines/plotting/colors/

const int NCMSColorScheme6 = 6;
const char* CMSColorScheme6Codes[NCMSColorScheme6] = {"#5790fc", "#f89c20", "#e42536", "#964a8b", "#9c9ca1", "#7a21dd"};

const int NCMSColorScheme10 = 10;
const char* CMSColorScheme10Codes[NCMSColorScheme10] = {"#3f90da", "#ffa90e", "#bd1f01", "#94a4a2", "#832db6", "#a96b59", "#e76300", "#b9ac70", "#717581", "#92dadd"};

const int CMSColorLightBlue = TColor::GetColor("#3f90da");
const int CMSColorOrange = TColor::GetColor("#ffa90e");
const int CMSColorRed = TColor::GetColor("#bd1f01");
const int CMSColorLightGray = TColor::GetColor("#94a4a2");
const int CMSColorPurple = TColor::GetColor("#832db6");
const int CMSColorBrown = TColor::GetColor("#a96b59");
const int CMSColorDarkGray = TColor::GetColor("#717581");

/// Colors picked up from the METbrewer colorblind-friendly palettes inspired by actual artworks at the MET in New-York https://github.com/BlakeRMills/MetBrewer

// "Sailing Boats Returning to Yabase", Lake Biwa, 1835, Utagawa Hiroshige https://github.com/BlakeRMills/MetBrewer?tab=readme-ov-file#hiroshige

const int HiroshigeLightRed = TColor::GetColor("#e76254");
const int HiroshigeOrange = TColor::GetColor("#ef8a47");
const int HiroshigeLightOrange = TColor::GetColor("#f7aa58");
const int HiroshigeYellow = TColor::GetColor("#ffd06f");
const int HiroshigeIceBlue = TColor::GetColor("#aadce0");
const int HiroshigeLightBlue = TColor::GetColor("#72bcd5");
const int HiroshigeBlue = TColor::GetColor("#528fad");
const int HiroshigeGrayBlue = TColor::GetColor("#376795");
const int HiroshigeNightBlue = TColor::GetColor("#1e466e");

// "Red and Yellow Cliffs", 1940, Georgia O'Keeffe https://github.com/BlakeRMills/MetBrewer?tab=readme-ov-file#okeeffe2
const int NOKeffeCliffsColors = 7;

const char* OKeeffeCliffsPalette[] = {"#fbe3c2", "#f2c88f", "#ecb27d", "#e69c6b", "#d37750", "#b9563f", "#92351e"};

// "Dragon Robe", 1998, Vivienne Tam https://github.com/BlakeRMills/MetBrewer/tree/main?tab=readme-ov-file#tam
const int NTamDragonColors = 8;

const char* TamDragonPalette[] = {"#ffd353", "#ffb242", "#ef8737", "#de4f33", "#bb292c", "#9f2d55", "#62205f", "#341648"};
const int TamDragonDarkPurple = TColor::GetColor("#341648");
const int TamDragonPurple = TColor::GetColor("#62205f");
const int TamDragonRedViolet = TColor::GetColor("#9f2d55");
const int TamDragonRed = TColor::GetColor("#bb292c");
const int TamDragonOrange = TColor::GetColor("#de4f33");
const int TamDragonLightOrange = TColor::GetColor("#ef8737");
const int TamDragonOrangeYellow = TColor::GetColor("#ffb242");
const int TamDragonYellow = TColor::GetColor("#ffd353");

const int TamDragonPaletteRed[] = {52, 98, 159, 197, 222, 239, 255, 255};
const int TamDragonPaletteGreen[] = {22, 32, 45, 41, 79, 135, 178, 211};
const int TamDragonPaletteBlue[] = {72, 95, 85, 44, 51, 55, 66, 83};

// Klimt's gold, inspired from "The Kiss"
const int NKlimtGoldColors = 5;

// RGB values
const int KlimtGoldPaletteRed[] = {128, 164, 201, 237, 255};
const int KlimtGoldPaletteGreen[] = {91, 126, 162, 197, 225};
const int KlimtGoldPaletteBlue[] = {16, 27, 39, 49, 105};

/// Palettes for multi-dimensional distributions

void SetColorPalette(std::string styleName = "cividis", bool invert = false, int nContours = 255) {
	// palettes defined in ROOT
	if (styleName == "cividis") {
		gStyle->SetPalette(kCividis);
	} else if (styleName == "viridis") {
		gStyle->SetPalette(kViridis);
	} else if (styleName == "sunset") {
		gStyle->SetPalette(kSunset);
	}

	// custom palettes based on the color schemes defined above

	else if (styleName == "CMS6") {
		const int nColors = NCMSColorScheme6;
		int colors[nColors];

		for (int i = 0; i < nColors; ++i) {
			colors[i] = TColor::GetColor(CMSColorScheme6Codes[i]);
		}
		gStyle->SetPalette(nColors, colors); // will not be a linear gradient though, too complicated to convert from hexadecimal to RBG codes
	}

	else if (styleName == "KlimtGold") {
		const int nColors = NKlimtGoldColors;

		double red[nColors], green[nColors], blue[nColors], stops[nColors];
		for (int i = 0; i < nColors; ++i) {
			red[i] = KlimtGoldPaletteRed[i] / 255.;
			green[i] = KlimtGoldPaletteGreen[i] / 255.;
			blue[i] = KlimtGoldPaletteBlue[i] / 255.;
			stops[i] = (double)i / (nColors - 1);
		}

		TColor::CreateGradientColorTable(nColors, stops, red, green, blue, nContours);
	}

	else if (styleName == "OKeeffeCliffs") {
		const int nColors = NOKeffeCliffsColors;
		int colors[nColors];

		for (int i = 0; i < nColors; ++i) {
			colors[i] = TColor::GetColor(OKeeffeCliffsPalette[i]);
		}
		gStyle->SetPalette(nColors, colors); // will not be a linear gradient though, too complicated to convert from hexadecimal to RBG codes
	}

	else if (styleName == "TamDragon") {
		const int nColors = NTamDragonColors;

		double red[nColors], green[nColors], blue[nColors], stops[nColors];
		for (int i = 0; i < nColors; ++i) {
			red[i] = TamDragonPaletteRed[i] / 255.;
			green[i] = TamDragonPaletteGreen[i] / 255.;
			blue[i] = TamDragonPaletteBlue[i] / 255.;
			stops[i] = (double)i / (nColors - 1);
		}

		TColor::CreateGradientColorTable(nColors, stops, red, green, blue, nContours);
	}

	else {
		std::cout << "Warning: the name provided does not correspond to any of the palette style name defined in ColorSchemes! Will set it to ROOT cividis by default...\n";
		gStyle->SetPalette(kCividis);
	}

	// gradient smoothing
	gStyle->SetNumberContours(nContours);

	if (invert) TColor::InvertPalette();
}
