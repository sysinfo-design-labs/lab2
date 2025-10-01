package com.example;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Main {

    private static final double TOL = 1e-9;
    private static final int MAX_BISECT_ITERS = 200;
    private static final double EPS = 1e-12;

    // Уравнение f(B) = 0 для поиска B
    private static double fOfB(double B, double[] Xi, double sumXi) {
        int n = Xi.length;
        double denom = 0.0;
        double sumInv = 0.0;

        for (int i = 1; i <= n; i++) {
            double term = B - i + 1.0;
            if (term <= 0) return Double.NaN;
            sumInv += 1.0 / term;
            denom += term * Xi[i - 1];
        }

        if (Math.abs(denom) < EPS) return Double.NaN;
        return sumInv - (n * sumXi) / denom;
    }

    // Поиск B методом бисекции
    private static double findB(double[] Xi) {
        int n = Xi.length;
        double sumXi = Arrays.stream(Xi).sum();

        double left = n + 1e-8;
        double fLeft = fOfB(left, Xi, sumXi);
        double right = left + 1.0;
        double fRight = fOfB(right, Xi, sumXi);

        int expansion = 0;
        while ((Double.isNaN(fLeft) || Double.isNaN(fRight) || fLeft * fRight > 0) && expansion < 1000) {
            right = left + (expansion + 1) * 2.0;
            fRight = fOfB(right, Xi, sumXi);
            expansion++;
            if (right > n + 1e6) break;
        }

        if (Double.isNaN(fLeft) || Double.isNaN(fRight) || fLeft * fRight > 0) {
            return n + 1.0; // fallback
        }

        double a = left, b = right;
        double fa = fOfB(a, Xi, sumXi), fb = fOfB(b, Xi, sumXi);
        for (int iter = 0; iter < MAX_BISECT_ITERS; iter++) {
            double m = 0.5 * (a + b);
            double fm = fOfB(m, Xi, sumXi);
            if (Double.isNaN(fm)) return m;
            if (Math.abs(fm) < TOL) return m;
            if (fa * fm <= 0) {
                b = m; fb = fm;
            } else {
                a = m; fa = fm;
            }
            if (Math.abs(b - a) < 1e-10) return 0.5 * (a + b);
        }
        return 0.5 * (a + b);
    }

    private static double computeK(double B, double[] Xi) {
        int n = Xi.length;
        double denom = 0.0;
        for (int i = 1; i <= n; i++) {
            denom += (B - i + 1.0) * Xi[i - 1];
        }
        if (Math.abs(denom) < EPS) return Double.NaN;
        return n / denom;
    }

    private static double meanTimeToNextFailure(double B, double K, int n) {
        if (B <= n) return Double.POSITIVE_INFINITY;
        double denom = K * (B - n);
        if (Math.abs(denom) < EPS) return Double.POSITIVE_INFINITY;
        return 1.0 / denom;
    }

    // исправленная формула для времени до окончания тестирования
    private static double timeToEndOfTesting(double B, double K, int n) {
        int remaining = (int)Math.floor(B) - n; // сколько ошибок осталось
        if (remaining <= 0) return 0.0;

        double sum = 0.0;
        for (int i = 1; i <= remaining; i++) {
            sum += 1.0 / i;
        }
        return sum / K;
    }

    // Парсинг чисел из строки
    private static double[] parseNumbersFromString(String s) {
        if (s == null) return new double[0];
        s = s.replace(',', '.');
        Pattern p = Pattern.compile("[+-]?\\d*\\.?\\d+");
        Matcher m = p.matcher(s);
        List<Double> list = new ArrayList<>();
        while (m.find()) {
            try {
                list.add(Double.parseDouble(m.group()));
            } catch (NumberFormatException ignored) {}
        }
        return list.stream().mapToDouble(Double::doubleValue).toArray();
    }

    // Главная функция
    public static void main(String[] args) {
        Scanner sc = new Scanner(System.in);
        try {
            System.out.println("Модель Джелинского–Моранды");
            System.out.println("Введите интервалы Xi (через пробел или запятую):");

            String line = sc.nextLine();
            double[] Xi = parseNumbersFromString(line);

            if (Xi.length == 0) {
                System.err.println("Ошибка: не найдено ни одного числа.");
                return;
            }

            int n = Xi.length;
            double sumXi = Arrays.stream(Xi).sum();

            System.out.printf("Число ошибок n = %d%n", n);
            System.out.printf("Сумма интервалов = %.6f%n", sumXi);

            double B = findB(Xi);
            System.out.printf("Оценка общего числа ошибок (B) = %.6f%n", B);

            double K = computeK(B, Xi);
            System.out.printf("Коэффициент пропорциональности (K) = %.8f%n", K);

            double Xnext = meanTimeToNextFailure(B, K, n);
            System.out.printf("Среднее время до следующей ошибки X(n+1) = %.6f  час%n", Xnext);

            double Tend = timeToEndOfTesting(B, K, n);
            System.out.printf("Время до окончания тестирования Tконец = %.6f час", Tend);

        } catch (Exception ex) {
            System.err.println("Ошибка: " + ex.getMessage());
        }
    }
}
