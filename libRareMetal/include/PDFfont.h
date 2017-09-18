////////////////////////////////////////////////////////////////////// 
// janpdf/PDFfont.h 
// (c) 2000-2007 Goncalo Abecasis (c) 2002-2007 Jan Wigginton
// 
// This file is distributed as part of the PEDSTATS source code package   
// and may not be redistributed in any form, without prior written    
// permission from the author. Permission is granted for you to       
// modify this file for your own personal use, but modified versions  
// must retain this copyright notice and must not be distributed.     
// 
// Permission is granted for you to use this file to compile PEDSTATS.    
// 
// All computer programs have bugs. Use this file at your own risk.   
// 
// Tuesday December 18, 2007
// 
 
#ifndef __PDFFONT_H__
#define __PDFFONT_H__

enum PDFFonts { fTimes = 0, fHelvetica = 1, fCourier = 2,
                fSymbol = 3, fZapfDingbats = 4 };

class PDF;

class PDFFont
   {
   private:
      static char * fontNames[14];
      static int  metrics[14][256];

      bool selectedFonts[14];
      int  fontDictionary;

      PDF & pdf;

   public:
      PDFFont(PDF & parent);

      // Functions for managing font dictionary
      void WriteResources();
      void WriteDictionary();
      void SelectFont(int fontid, double size);

      // Returns ID 0-14 for selected font
      int    GetFontID  (PDFFonts font, bool bold, bool italic );

      // Get Type 1 font name for a specific font
      char * GetFontName(PDFFonts font, bool bold, bool italic );
      char * GetFontName(int fontid);

      // Calculate width for a text string
      int  TextWidth(int font, const char * string);

      // Marks font for inclusion in font dictionary
      void MarkFont(int fontid);

   };

#endif


 
 
 
 
 
 
 
 
 
 
