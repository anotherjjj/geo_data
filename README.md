# Methylation Dataset Creator

Скрипт для сбора данных о метилировании ДНК из базы данных GEO.

## Описание

Создает три CSV файла:
- `methylation_gse_ids.csv` - список идентификаторов datasets
- `methylation_datasets_info.csv` - информация о datasets (название, описание, организм, платформа)
- `methylation_samples_metadata.csv` - метаданные образцов (возраст, пол, ткань)

## Использование

```bash
# Создать датасет с параметрами по умолчанию
python tablecreator.py

```
