from openpyxl import Workbook
from openpyxl.drawing.image import Image
from openpyxl.styles import PatternFill


def make_exel_file(params, file_name, best_correlations):
    column_names = ['pKa', 'dih angle', 'fitness']
    correlation_column_names = ['increase pKa', 'decrease pKa', 'increase dih angle', 'dec dih angle']
    wb = Workbook()
    ws = wb.active
    letter_first = 65
    step = 8    # 4 - это отступ в один столбец между
    # Заголовки колонок
    for column_idx in range(len(params)):
        ws[f"{chr(letter_first)}1"] = f"Parameter {column_names[column_idx]}"
        ws[f"{chr(letter_first)}2"] = "mutations"
        ws[f"{chr(letter_first + 1)}2"] = "the first par"
        ws[f"{chr(letter_first + 2)}2"] = "the second par"

        ws[f"{chr(letter_first + 3)}2"] = correlation_column_names[0]
        ws[f"{chr(letter_first + 4)}2"] = correlation_column_names[1]
        ws[f"{chr(letter_first + 5)}2"] = correlation_column_names[2]
        ws[f"{chr(letter_first + 6)}2"] = correlation_column_names[3]
        letter_first += step
    colors = ["7cedac", "97e6ed", "dd76ef", "9f88ec"]
    for i, slices in enumerate(params):
        for j, muts in enumerate(slices):
            cell = ws.cell(row=j + 3, column=i*step + 1)
            cell.value = muts[0]
            cell = ws.cell(row=j + 3, column=i * step + 2)
            cell.value = round(muts[1][0], 2)
            cell = ws.cell(row=j + 3, column=i * step + 3)
            cell.value = round(muts[1][1], 2)
            for col_idx, cor in enumerate(best_correlations):
                if any([x[0] in muts[0] for x in cor]):
                    fill = PatternFill(start_color=colors[col_idx], end_color=colors[col_idx], fill_type="solid")
                    cell = ws.cell(row=j + 3, column=i * step + 4 + col_idx)
                    cell.fill = fill


    for idx, coord in zip([1, 2, 3], [f'Y2', 'Y20', 'AK2']):
        image_path = f'his_{idx}.png'
        img = Image(image_path)
        ws.add_image(img, coord)
        # Сохранение изменений
        wb.save(f"output_{file_name}.xlsx")