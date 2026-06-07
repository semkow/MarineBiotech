# ==============================================================================
# SKRYPT DO BIOMARKERÓW I METABOLOMIKI: SPAJANIE DUPLIKATÓW POD GNPS2 / AAMN
# ==============================================================================

# 1. ŁADOWANIE BIBLIOTEK
# Jeśli nie masz pakietu tidyverse, najpierw odkomentuj i uruchom linijkę poniżej:
# install.packages("tidyverse")
library(tidyverse)

# 2. OKIENKO WYBORU PLIKU WEJŚCIOWEGO
# Otworzy się okno systemowe, w którym wskazujesz oczyszczony plik .csv
sciezka_pliku <- file.choose()
cat("Pomyślnie załadowano plik:", sciezka_pliku, "\n")

# Wczytanie danych do pamięci R
tabela_wejsciowa <- read.csv(sciezka_pliku, check.names = FALSE)

# 3. DEFINICJA FUNKCJI ŁĄCZĄCEJ WIERSZE
spajaj_duplikaty <- function(dane) {
  
  # Sortowanie danych chromatograficznych dla porządku w pętli
  dane <- dane %>% arrange(`row m/z`, `row retention time`)
  
  # Automatyczne wykrywanie kolumn z próbkami (szukamy frazy "Peak area")
  kolumny_peak_area <- grep("Peak area", names(dane), value = TRUE)
  
  # Stworzenie wektora indeksów do grupowania
  numery_grup <- 1:nrow(dane)
  
  # Pętla sprawdzająca Twoje precyzyjne kryteria odległości masowej i czasowej
  for (i in 1:(nrow(dane) - 1)) {
    for (j in (i + 1):nrow(dane)) {
      
      roznica_mz <- abs(dane$`row m/z`[i] - dane$`row m/z`[j])
      roznica_rt <- abs(dane$`row retention time`[i] - dane$`row retention time`[j])
      
      # KRYTERIUM 1: Różnica m/z <= 0.01 ORAZ różnica RT <= 0.02
      warunek_bliski <- (roznica_mz <= 0.5) && (roznica_rt <= 0.1)
      
      # KRYTERIUM 2: Różnica m/z = 1 (+/- 0.01) ORAZ różnica RT <= 0.1
      warunek_izotop_addukt <- (abs(roznica_mz - 1) <= 0.1) && (roznica_rt <= 0.1)
      
      # Jeśli którykolwiek z Twoich dwóch warunków jest spełniony -> spajamy wiersze
      if (warunek_bliski || warunek_izotop_addukt) {
        najnizszy_indeks <- min(numery_grup[i], numery_grup[j])
        numery_grup[numery_grup == numery_grup[i]] <- najnizszy_indeks
        numery_grup[numery_grup == numery_grup[j]] <- najnizszy_indeks
      }
    }
  }
  
  # Przypisujemy wyliczone grupy do tabeli
  dane$identyfikator_grupy <- numery_grup
  
  # Główna agregacja matematyczna danych
  tabela_wynikowa <- dane %>%
    group_by(identyfikator_grupy) %>%
    summarise(
      `row ID` = min(`row ID`),                                # Wybiera mniejsze ID spośród wierszy
      `row m/z` = mean(`row m/z`),                             # Wylicza średnią masę m/z
      `row retention time` = mean(`row retention time`),       # Wylicza średni czas retencji RT
      
      # Sumowanie wartości wyłącznie w kolumnach z intensywnościami "Peak area"
      across(all_of(kolumny_peak_area), ~ sum(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    # Usunięcie technicznej kolumny grupującej przed eksportem
    select(-identyfikator_grupy)
  
  return(tabela_wynikowa)
}

# 4. URUCHOMIENIE ALGORYTMU
cat("Trwa przetwarzanie danych i scalanie duplikatów... Proszę czekać.\n")
tabela_czysta <- spajaj_duplikaty(tabela_wejsciowa)

# 5. AUTOMATYCZNY ZAPIS PLIKU W TYM SAMYM FOLDERZE
folder_docelowy <- dirname(sciezka_pliku)
nazwa_nowego_pliku <- paste0("CLEAN_", basename(sciezka_pliku))
pelna_sciezka_zapisu <- file.path(folder_docelowy, nazwa_nowego_pliku)

write.csv(tabela_czysta, pelna_sciezka_zapisu, row.names = FALSE)

# KOMUNIKAT KOŃCOWY
cat("\n==================================================================\n")
cat("SUKCES! Duplikaty zostały pomyślnie sklejone.\n")
cat("Nowy plik został zapisany jako:", nazwa_nowego_pliku, "\n")
cat("Znajdziesz go w tym samym folderze, co plik źródłowy.\n")
cat("==================================================================\n")