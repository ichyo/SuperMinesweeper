package com.topcoder.marathon;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.UncheckedIOException;

class ErrorReader {
    private final BufferedReader errorStream;
    private final BufferedWriter errorWriter;
    private final StringBuilder sb = new StringBuilder();
    private final boolean printMessages;
    private static final int maxLength = 10_000_000;

    public ErrorReader(BufferedReader errorStream, boolean printMessages, BufferedWriter errorWriter) {
        this.errorStream = errorStream;
        this.printMessages = printMessages;
        this.errorWriter = errorWriter;
    }

    public void readAndWrite() throws IOException {
        int ch;
        final StringBuilder buffer = new StringBuilder();
        while (errorStream.ready() && (ch = errorStream.read()) >= 0) {
            buffer.append((char)ch);
        }
        if (buffer.length() > 0) {
            write(buffer.toString());
        }
    }

    private void write(String s) throws IOException {
        if (sb.length() < maxLength) sb.append(s);
        if (printMessages) {
            System.out.print(s);
            System.out.flush();
        }
        if (errorWriter != null) {
            errorWriter.write(s);
            errorWriter.flush();
        }
    }

    public String getOutput() {
        return sb.toString();
    }

    public void close() {
        try {
            readAndWrite();
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }

        try {
            if (errorStream != null) errorStream.close();
        } catch (Exception e) {
        }
        try {
            if (errorWriter != null) errorWriter.close();
        } catch (Exception e) {
        }
    }
}