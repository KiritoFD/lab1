// Include necessary headers
// ...

// Renamed the main function to avoid conflicts
int generator_main(int argc, char* argv[]) {
    // Original main function contents
    // ...
    return 0;
}

#ifdef GENERATOR_BUILD
// Only include this main function when specifically building the generator
int main(int argc, char* argv[]) {
    return generator_main(argc, argv);
}
#endif
