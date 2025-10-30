# ==============================================================================
# Load IP-MS results file (sent as a .xlsx).
# ==============================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

load_ip_ms_data <- function(
    filename,
    sheet
) {
  dat <- readxl::read_excel(
    path = filename,
    sheet = sheet
  )
  return(dat)
}
